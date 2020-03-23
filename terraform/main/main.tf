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
  name = "allow_internal"
  ingress {
    from_port = 0
    to_port   = 0
    protocol  = "-1"
    self      = true
  }
}

resource "aws_s3_bucket" "work" {
  bucket_prefix = "tf-serratus-work-"
  force_destroy = true
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

module "splitter" {
  source = "../worker"

  up                 = var.up
  dev_cidrs          = var.dev_cidrs
  security_group_ids = [aws_security_group.internal.id]
  instance_type      = "t3.small"
  spot_price         = 0.007
  volume_size        = 50
  s3_bucket          = aws_s3_bucket.work.bucket
  s3_prefix          = "fq-blocks"
  dockerhub_account  = var.dockerhub_account
  image_name         = "serratus-dl"
  key_name           = var.key_name
  scheduler          = "${module.scheduler.public_dns}:${var.scheduler_port}"
  options            = "-k s3://${aws_s3_bucket.work.bucket}/fq-blocks"
}

output "scheduler_dns" {
  value = module.scheduler.public_dns
}
