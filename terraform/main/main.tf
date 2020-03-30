variable "aws_region" {
  type    = string
  default = "us-east-1"
}

variable "up" {
  type    = bool
  default = true
}

variable "dev_cidrs" {
  type = set(string)
}

variable "key_name" {
  description = "Name of the previously created EC2 key pair, for SSH access"
  type        = string
}

variable "dockerhub_account" {
  type = string
}

variable "scheduler_port" {
  type  = number
  default = 8000
}

provider "aws" {
  region = var.aws_region
}

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

resource "aws_s3_bucket" "work" {
  bucket_prefix = "tf-serratus-work-"
  force_destroy = true

  tags = {
    "project": "serratus"
    "component": "serratus-scheduler"
  }
}

module "scheduler" {
  source = "../scheduler"

  dev_cidrs          = var.dev_cidrs
  security_group_ids = [aws_security_group.internal.id]
  instance_type      = "t3.nano"
  dockerhub_account  = var.dockerhub_account
  key_name           = var.key_name
  scheduler_port     = var.scheduler_port
}

module "download" {
  source = "../worker"

  up                 = var.up
  dev_cidrs          = var.dev_cidrs
  security_group_ids = [aws_security_group.internal.id]
  instance_type      = "c5.large"
  spot_price         = 0.04
  s3_bucket          = aws_s3_bucket.work.bucket
  s3_prefix          = "fq-blocks"
  dockerhub_account  = var.dockerhub_account
  image_name         = "serratus-dl"
  workers            = 2
  key_name           = var.key_name
  scheduler          = "${module.scheduler.public_dns}:${var.scheduler_port}"
  options            = "-k ${aws_s3_bucket.work.bucket}"
}

module "align" {
  source = "../worker"

  up                 = var.up
  dev_cidrs          = var.dev_cidrs
  security_group_ids = [aws_security_group.internal.id]
  instance_type      = "c5.large" # c5.large
  spot_price         = 0.04
  s3_bucket          = aws_s3_bucket.work.bucket
  s3_delete_prefix   = "fq-blocks"
  s3_prefix          = "bam-blocks"
  dockerhub_account  = var.dockerhub_account
  image_name         = "serratus-align"
  key_name           = var.key_name
  scheduler          = "${module.scheduler.public_dns}:${var.scheduler_port}"
  options            = "-k ${aws_s3_bucket.work.bucket}"
}

module "merge" {
  source = "../worker"

  up                 = var.up
  dev_cidrs          = var.dev_cidrs
  security_group_ids = [aws_security_group.internal.id]
  instance_type      = "t3.small"
  spot_price         = 0.007
  s3_bucket          = aws_s3_bucket.work.bucket
  s3_delete_prefix   = "bam-blocks"
  s3_prefix          = "out"
  dockerhub_account  = var.dockerhub_account
  image_name         = "serratus-merge"
  key_name           = var.key_name
  scheduler          = "${module.scheduler.public_dns}:${var.scheduler_port}"
  options            = "-k ${aws_s3_bucket.work.bucket} -b s3://${aws_s3_bucket.work.bucket}/out"
}

output "scheduler_dns" {
  value = module.scheduler.public_dns
}
