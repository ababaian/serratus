provider "aws" {
  region = "us-east-1"
}

variable "up" {
  type    = bool
  default = true
}

variable "dev_cidrs" {
  type    = set(string)
  default = ["173.181.15.58/32"] # Jeff's house
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
}

module "splitter" {
  source = "../worker"

  up                 = var.up
  dev_cidrs          = var.dev_cidrs
  scheduler_dns      = module.scheduler.public_dns
  security_group_ids = [aws_security_group.internal.id]
  instance_type      = "t3.small"
  spot_price         = 0.007
  volume_size        = 50
  s3_bucket          = aws_s3_bucket.work.bucket
  s3_prefix          = "fq-blocks"
}

output "scheduler_dns" {
  value = module.scheduler.public_dns
}
