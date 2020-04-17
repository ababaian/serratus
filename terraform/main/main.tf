///////////////////////////////////////////////////////////
// SERRATUS MAIN TERRAFORM CONFIG
///////////////////////////////////////////////////////////
// VARIABLES ##############################

variable "aws_region" {
  type    = string
  default = "us-east-1"
}

variable "dl_size" {
  type    = number
  default = 0
  description = "Default number of downloader nodes (ASG)"
}

variable "align_size" {
  type    = number
  default = 0
  description = "Default number of aligner nodes (ASG)"
}

variable "dev_cidrs" {
  type = set(string)
}

variable "key_name" {
  description = "Name of the previously created EC2 key pair, for SSH access"
  type        = string
}

variable "dockerhub_account" {
  type        = string
  description = "Dockerhub account name, where images will be pulled from"
}

variable "scheduler_port" {
  type  = number
  default = 8000
}

// PROVIDER/AWS ##############################
provider "aws" {
  version     = "~> 2.49"
  region      = var.aws_region
}

provider "local" {}

resource "aws_security_group" "internal" {
  name = "serratus-internal"
  ingress {
    # Allow all internal traffic between workers
    from_port = 0
    to_port   = 0
    protocol  = "-1"
    self      = true
  }
  ingress {
    # Allow dev SSH access from a minimal IP range
    from_port   = 22
    to_port     = 22
    protocol    = "tcp"
    cidr_blocks = var.dev_cidrs
  }
  egress {
    # Required to download docker images and yum-update.
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

// MODULES ##############################

// Working S3 storage for Serratus
module "work_bucket" {
  source   = "../bucket"

  prefixes = ["fq-blocks", "bam-blocks", "out"]
}

// Cluster scheduler and task manager
module "scheduler" {
  source             = "../scheduler"
  
  security_group_ids = [aws_security_group.internal.id]
  key_name           = var.key_name
  instance_type      = "t3.nano"
  dockerhub_account  = var.dockerhub_account
  scheduler_port     = var.scheduler_port
}

// Cluster monitor
module "monitoring" {
  source             = "../monitoring"

  security_group_ids = [aws_security_group.internal.id]
  key_name           = var.key_name
  scheduler_ip       = module.scheduler.private_ip
  instance_type      = "t3.medium"
}

// Serratus-dl
module "download" {
  source             = "../worker"

  desired_size       = 0
  max_size           = 256
  dev_cidrs          = var.dev_cidrs
  security_group_ids = [aws_security_group.internal.id]
  instance_type      = "c5.large"
  spot_price         = 0.04
  s3_bucket          = module.work_bucket.name
  s3_prefix          = "fq-blocks"
  dockerhub_account  = var.dockerhub_account
  image_name         = "serratus-dl"
  workers            = 2
  key_name           = var.key_name
  scheduler          = "${module.scheduler.public_dns}:${var.scheduler_port}"
  options            = "-k ${module.work_bucket.name}"
}

// Serratus-align
module "align" {
  source             = "../worker"

  desired_size       = 0
  max_size           = 256
  dev_cidrs          = var.dev_cidrs
  security_group_ids = [aws_security_group.internal.id]
  instance_type      = "c5.large" # c5.large
  spot_price         = 0.04
  s3_bucket          = module.work_bucket.name
  s3_delete_prefix   = "fq-blocks"
  s3_prefix          = "bam-blocks"
  dockerhub_account  = var.dockerhub_account
  image_name         = "serratus-align"
  key_name           = var.key_name
  scheduler          = "${module.scheduler.public_dns}:${var.scheduler_port}"
  options            = "-k ${module.work_bucket.name} -a bowtie2"
}

//Serratus-merge
module "merge" {
  source             = "../worker"

  desired_size       = 1
  max_size           = 10
  dev_cidrs          = var.dev_cidrs
  security_group_ids = [aws_security_group.internal.id]
  instance_type      = "t3.small"
  spot_price         = 0.007
  s3_bucket          = module.work_bucket.name
  s3_delete_prefix   = "bam-blocks"
  s3_prefix          = "out"
  dockerhub_account  = var.dockerhub_account
  image_name         = "serratus-merge"
  key_name           = var.key_name
  scheduler          = "${module.scheduler.public_dns}:${var.scheduler_port}"
  options            = "-k ${module.work_bucket.name} -b s3://${module.work_bucket.name}/out"
}

// RESOURCES ##############################
// Controller scripts created locally

resource "local_file" "hosts" {
  filename = "${path.module}/serratus-hosts"
  file_permission = 0666
  content = <<-EOF
    aws_monitor ${module.monitoring.public_dns}
    aws_scheduler ${module.scheduler.public_dns}
  EOF
}

resource "local_file" "create_tunnel" {
  filename = "${path.module}/create_tunnels.sh"
  file_permission = 0777
  content = <<-EOF
    #!/bin/bash
    set -eu
    echo "Tunnels created:"
    ssh -Nf -L 3000:localhost:3000 ec2-user@${module.monitoring.public_dns}
    echo "    localhost:3000 -- grafana"
    ssh -Nf -L 9090:localhost:9090 ec2-user@${module.monitoring.public_dns}
    echo "    localhost:9090 -- prometheus"
    ssh -Nf -L 8000:localhost:8000 ec2-user@${module.scheduler.public_dns}
    echo "    localhost:8000 -- scheduler"
  EOF
}

resource "local_file" "dl_set_capacity" {
  filename = "${path.module}/dl_set_capacity.sh"
  file_permission = 0777
  content = <<-EOF
    #!/bin/bash
    set -eux
    aws autoscaling set-desired-capacity \
      --auto-scaling-group-name ${module.download.asg_name} \
      --desired-capacity $1
  EOF
}

resource "local_file" "align_set_capacity" {
  filename = "${path.module}/align_set_capacity.sh"
  file_permission = 0777
  content = <<-EOF
    #!/bin/bash
    set -eux
    aws autoscaling set-desired-capacity \
      --auto-scaling-group-name ${module.align.asg_name} \
      --desired-capacity $1
  EOF
}

resource "local_file" "merge_set_capacity" {
  filename = "${path.module}/merge_set_capacity.sh"
  file_permission = 0777
  content = <<-EOF
    #!/bin/bash
    set -eux
    aws autoscaling set-desired-capacity \
      --auto-scaling-group-name ${module.merge.asg_name} \
      --desired-capacity $1
  EOF
}

// OUTPUT ##############################
output "help" {
  value = <<-EOF
    Run ${local_file.create_tunnel.filename} to create SSH tunnels for all services.
  EOF
}

output "scheduler_dns" {
  value = module.scheduler.public_dns
}

output "monitor_dns" {
  value = module.monitoring.public_dns
}

output "dl_asg_name" {
  value = module.download.asg_name
}
output "merge_asg_name" {
  value = module.merge.asg_name
}
output "align_asg_name" {
  value = module.align.asg_name
}
