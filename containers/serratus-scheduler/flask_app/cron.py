import os
import time
import math
from multiprocessing import Process

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


def set_asg_size(autoscaling, asg_name, constant, max_, num_jobs, name):
    desired_size = min(max_, math.ceil(constant * num_jobs))
    print("Adjusting {} ASG to {}".format(name, desired_size))
    try:
        autoscaling.set_desired_capacity(
            AutoScalingGroupName=asg_name,
            DesiredCapacity=desired_size,
        )
    except (autoscaling.exceptions.ScalingActivityInProgressFault,
            autoscaling.exceptions.ResourceContentionFault) as e:
        print("Error setting ASG size, {}".format(e))


def adjust_autoscaling_loop(app):
    autoscaling = boto3.client('autoscaling', region_name=app.config["AWS_REGION"])
    dl_asg = get_asg_name(autoscaling, "serratus-dl")
    align_asg = get_asg_name(autoscaling, "serratus-align")
    merge_asg = get_asg_name(autoscaling, "serratus-merge")
    while True:
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

        if config["DL_SCALING_ENABLE"]:
            constant = float(config["DL_SCALING_CONSTANT"])
            max_ = int(config["DL_SCALING_MAX"])
            set_asg_size(
                autoscaling, dl_asg, constant, max_, num_dl_jobs, "dl"
            )
        if config["ALIGN_SCALING_ENABLE"]:
            constant = float(config["ALIGN_SCALING_CONSTANT"])
            max_ = int(config["ALIGN_SCALING_MAX"])
            set_asg_size(
                autoscaling, align_asg, constant, max_, num_align_jobs, "align"
            )
        if config["MERGE_SCALING_ENABLE"]:
            constant = float(config["MERGE_SCALING_CONSTANT"])
            max_ = int(config["MERGE_SCALING_MAX"])
            set_asg_size(
                autoscaling, merge_asg, constant, max_, num_merge_jobs, "merge"
            )

        scale_interval = int(config["SCALING_INTERVAL"])
        print("ajust_autoscaling() finished.  Running again in {} seconds".format(scale_interval))
        time.sleep(scale_interval)


def clean_terminated_jobs_loop(app):
    while True:
        print("clear_terminated_jobs() triggered by interval")
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
        Process(target=clean_terminated_jobs_loop, args=(app,)).start()
        Process(target=adjust_autoscaling_loop, args=(app,)).start()
