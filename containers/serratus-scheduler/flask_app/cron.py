"""Cron system, runs periodic database cleanup jobs.

Run with `flask-3 cron`
"""
import time
import math
from threading import Thread

import boto3
import click
from flask import current_app
from flask.cli import with_appcontext
from prometheus_client import Gauge, start_http_server, CollectorRegistry
from sqlalchemy import or_

from . import db

CRON_REGISTRY = CollectorRegistry()

# Usage note: Prometheus requires special treatment if run with multiprocessing.
asg_desired_size_gauge = Gauge(
    "asg_desired_size", "Desired size of ASGs", ["asg"], registry=CRON_REGISTRY
)
asg_virt_size_gauge = Gauge(
    "asg_virt_size", "Desired size of ASGs", ["asg"], registry=CRON_REGISTRY
)


def get_asg_name(autoscaling, pattern):
    """Look through all available ASGs to find one matching a pattern"""
    paginator = autoscaling.get_paginator("describe_auto_scaling_groups").paginate()
    for page in paginator:
        for asg in page["AutoScalingGroups"]:
            name = asg["AutoScalingGroupName"]
            if pattern in name:
                return name
    raise RuntimeError('No ASG named "{}"'.format(pattern))


def set_asg_size(
    autoscaling,
    asg_name,
    constant,
    max_,
    num_jobs,
    name,
    virt_max_increase,
    _virt_sizes={},
):
    true_desired_size = min(max_, math.ceil(constant * num_jobs))

    # XXX: Creating too many instances causes problems in serratus-dl when
    # it tries to query SRA, maybe due to too many simultaneous queries?
    # We're fixing this by rate-limiting the increases to ASG size.
    # A better solution would be to implement exponential backoff where the
    # API is called, but that's buried somewhere in the code for fastq-dump.
    virt_size_max = _virt_sizes.get(name, 0) + virt_max_increase
    asg_desired_size = min(true_desired_size, virt_size_max)
    _virt_sizes[name] = asg_desired_size

    virt_asg = asg_desired_size != true_desired_size

    asg_desired_size_gauge.labels(name).set(true_desired_size)
    asg_virt_size_gauge.labels(name).set(asg_desired_size)

    try:
        autoscaling.set_desired_capacity(
            AutoScalingGroupName=asg_name, DesiredCapacity=asg_desired_size,
        )
        return virt_asg
    except (
        autoscaling.exceptions.ScalingActivityInProgressFault,
        autoscaling.exceptions.ResourceContentionFault,
        autoscaling.exceptions.ClientError,
    ) as e:
        print("Error setting ASG size, {}".format(e))
        return False


def adjust_autoscaling_loop(app):
    time.sleep(10)  # Give postgres a few seconds to start
    autoscaling = boto3.session.Session().client(
        "autoscaling", region_name=app.config["AWS_REGION"]
    )
    dl_asg = get_asg_name(autoscaling, "serratus-dl")
    align_asg = get_asg_name(autoscaling, "serratus-align")
    merge_asg = get_asg_name(autoscaling, "serratus-merge")

    while True:
        virt_asg = False
        with app.app_context():
            config = dict(db.get_config())

            session = db.get_session()
            num_dl_jobs = (
                session.query(db.Accession)
                .filter(
                    or_(db.Accession.state == "new", db.Accession.state == "splitting")
                )
                .count()
            )

            num_align_jobs = (
                session.query(db.Block)
                .filter(or_(db.Block.state == "new", db.Block.state == "aligning"))
                .count()
            )

            num_merge_jobs = (
                session.query(db.Accession)
                .filter(
                    or_(
                        db.Accession.state == "merge_wait",
                        db.Accession.state == "merging",
                    )
                )
                .count()
            )

        virt_max = config["VIRTUAL_ASG_MAX_INCREASE"]
        if config["DL_SCALING_ENABLE"]:
            constant = float(config["DL_SCALING_CONSTANT"])
            max_ = int(config["DL_SCALING_MAX"])
            virt_asg = virt_asg or set_asg_size(
                autoscaling, dl_asg, constant, max_, num_dl_jobs, "dl", virt_max,
            )
        if config["ALIGN_SCALING_ENABLE"]:
            constant = float(config["ALIGN_SCALING_CONSTANT"])
            max_ = int(config["ALIGN_SCALING_MAX"])
            virt_asg = virt_asg or set_asg_size(
                autoscaling,
                align_asg,
                constant,
                max_,
                num_align_jobs,
                "align",
                virt_max,
            )
        if config["MERGE_SCALING_ENABLE"]:
            constant = float(config["MERGE_SCALING_CONSTANT"])
            max_ = int(config["MERGE_SCALING_MAX"])
            virt_asg = virt_asg or set_asg_size(
                autoscaling,
                merge_asg,
                constant,
                max_,
                num_merge_jobs,
                "merge",
                virt_max,
            )

        scale_interval = int(config["SCALING_INTERVAL"])
        virt_interval = int(config["VIRTUAL_SCALING_INTERVAL"])
        if virt_asg:
            print(
                "virtual autoscaling finished.  Running again in {} seconds".format(
                    virt_interval
                )
            )
            time.sleep(virt_interval)
        else:
            print(
                "ajust_autoscaling() finished.  Running again in {} seconds".format(
                    scale_interval
                )
            )
            time.sleep(scale_interval)


def get_running_instances():
    ec2 = boto3.session.Session().client(
        "ec2", region_name=current_app.config["AWS_REGION"]
    )
    """Get a list of all EC2 instance IDs currently running"""
    for r in ec2.describe_instances()["Reservations"]:
        for instance in r["Instances"]:
            if instance["State"]["Name"] == "running":
                yield instance["InstanceId"]


def worker_to_instance_id(worker_id):
    # Example worker_id: "i-0a9a0d75781577180-7"
    # Desired instance id: "i-0a9a0d75781577180"
    return "-".join(worker_id.split("-")[:-1])


def check_and_clear(instances, table, active_state, new_state, name):
    """Check and reset jobs of a given type."""
    session = db.get_session()
    missing_instances = list()

    accessions = (
        session.query(table).filter(table.state == active_state).with_for_update().all()
    )

    count = 0
    for accession in accessions:
        if name == "dl":
            worker_id = accession.split_worker
        elif name == "merge":
            worker_id = accession.merge_worker
        elif name == "align":
            worker_id = accession.align_worker
        else:
            raise AssertionError("Invalid job type {}".format(name))

        instance_id = worker_to_instance_id(worker_id)

        if instance_id not in instances:
            accession.state = new_state
            missing_instances.append(instance_id)
            count += 1

    if missing_instances:
        print(
            "Reset jobs on {} {} instances, which were terminated:".format(
                len(missing_instances), name
            )
        )
        for instance in missing_instances:
            print("   {}".format(instance))

    if count:
        session.commit()

    return count


def clear_terminated_jobs():
    """Reset all jobs (dl, align, merge) which is in the running state but
    where the instance no longer exists.

    This should run inside of a session context, since the current DB doesn't
    handle transactions well.  What we should do is implement a global DB lock
    but I would need to test how that impacts performance."""
    instances = set(get_running_instances())

    check_and_clear(instances, db.Accession, "splitting", "new", "dl")
    check_and_clear(instances, db.Accession, "merging", "merge_wait", "merge")
    check_and_clear(instances, db.Block, "aligning", "new", "align")


def clean_terminated_jobs_loop(app):
    time.sleep(10)  # Give postgres a few seconds to start
    while True:
        with app.app_context():
            clear_interval = int(db.get_config_val("CLEAR_INTERVAL"))
            clear_terminated_jobs()

        print(
            "clear_terminated_jobs() finished.  Running again in {} seconds".format(
                clear_interval
            )
        )
        time.sleep(clear_interval)


@click.command("cron")
@with_appcontext
def cron():
    print("Creating background processes")
    app = current_app._get_current_object()
    start_http_server(9101, registry=CRON_REGISTRY)
    Thread(target=clean_terminated_jobs_loop, args=(app,)).start()
    Thread(target=adjust_autoscaling_loop, args=(app,)).start()


def register(app):
    app.cli.add_command(cron)
