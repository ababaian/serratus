import os
import time
import math
from threading import Thread

from prometheus_client import Gauge
from sqlalchemy import or_
import boto3

from . import db, jobs

def get_asg_name(autoscaling, pattern):
    """Look through all available ASGs to find one matching a pattern"""
    paginator = autoscaling.get_paginator("describe_auto_scaling_groups").paginate()
    for page in paginator:
        for asg in page["AutoScalingGroups"]:
            name = asg["AutoScalingGroupName"]
            if pattern in name:
                return name
    raise RuntimeError('No ASG named "{}"'.format(pattern))

# Usage note: Prometheus requires special treatment if run with multiprocessing.
asg_desired_size_gauge = Gauge('asg_desired_size', 'Desired size of ASGs', ['asg'])
asg_virt_size_gauge = Gauge('asg_virt_size', 'Desired size of ASGs', ['asg'])

def set_asg_size(
        autoscaling, asg_name, constant, max_, num_jobs, name,
        virt_max_increase, _virt_sizes={}):
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
            AutoScalingGroupName=asg_name,
            DesiredCapacity=asg_desired_size,
        )
        return virt_asg
    except (autoscaling.exceptions.ScalingActivityInProgressFault,
            autoscaling.exceptions.ResourceContentionFault,
            autoscaling.exceptions.ClientError) as e:
        print("Error setting ASG size, {}".format(e))
        return False



def adjust_autoscaling_loop(app):
    autoscaling = boto3.client('autoscaling', region_name=app.config["AWS_REGION"])
    dl_asg = get_asg_name(autoscaling, "serratus-dl")
    align_asg = get_asg_name(autoscaling, "serratus-align")
    merge_asg = get_asg_name(autoscaling, "serratus-merge")

    while True:
        virt_asg = False
        with app.app_context():
            config = dict(db.get_config())

            session = db.get_session()
            num_dl_jobs = session.query(db.Accession)\
                .filter(or_(db.Accession.state == "new",
                            db.Accession.state == "splitting"))\
                .count()

            num_align_jobs = session.query(db.Block)\
                .filter(or_(db.Block.state == "new",
                            db.Block.state == "aligning"))\
                .count()

            num_merge_jobs = session.query(db.Accession)\
                .filter(or_(db.Accession.state == "merge_wait",
                            db.Accession.state == "merging"))\
                .count()

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
                autoscaling, align_asg, constant, max_, num_align_jobs, "align", virt_max,
            )
        if config["MERGE_SCALING_ENABLE"]:
            constant = float(config["MERGE_SCALING_CONSTANT"])
            max_ = int(config["MERGE_SCALING_MAX"])
            virt_asg = virt_asg or set_asg_size(
                autoscaling, merge_asg, constant, max_, num_merge_jobs, "merge", virt_max,
            )

        scale_interval = int(config["SCALING_INTERVAL"])
        virt_interval = int(config["VIRTUAL_SCALING_INTERVAL"])
        if virt_asg:
            print("virtual autoscaling finished.  Running again in {} seconds".format(virt_interval))
            time.sleep(virt_interval)
        else:
            print("ajust_autoscaling() finished.  Running again in {} seconds".format(scale_interval))
            time.sleep(scale_interval)


def clean_terminated_jobs_loop(app):
    while True:
        with app.app_context():
            jobs.clear_terminated_jobs()
            clear_interval = int(db.get_config_val("CLEAR_INTERVAL"))
        print("clear_terminated_jobs() finished.  Running again in {} seconds"
              .format(clear_interval))
        time.sleep(clear_interval)


def register(app):
    # Prevent this scheduler from running two instances in debug mode.
    # See https://stackoverflow.com/a/25519547 for details
    if not app.debug or os.environ.get('WERKZEUG_RUN_MAIN') == 'true':
        print("Creating new process")
        Thread(target=clean_terminated_jobs_loop, args=(app,), daemon=True).start()
        Thread(target=adjust_autoscaling_loop, args=(app,), daemon=True).start()
