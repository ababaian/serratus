variable "up" {
  description = "Spin up actual instances"
  type        = bool
  default     = true
}

variable "scheduler_port" {
  description = "The port the worker will use to contact the scheduler"
  type        = number
  default     = 8000
}

variable "scheduler_dns" {
  description = "IP Address of the scheduler"
  type        = string
}

variable "instance_type" {
  description = "Type of node to use for the workers"
  type        = string
}

variable "spot_price" {
  type    = number
  default = 0.0016
}

variable "volume_size" {
  type    = number
}

variable "allow_ssh" {
  description = "Allow SSH access to the nodes"
  type        = bool
  default     = true
}

variable "dev_cidrs" {
  description = "Remote IP Address, for SSH, HTTP, etc access"
  type        = set(string)
}

variable "image_name" {
   description = "Docker image to run once the container is started"
   type        = string
   default     = "serratus-dl"
}

variable "s3_bucket" {
   type        = string
   description = "Name of the S3 bucket to store temporary data"
}

variable "s3_prefix" {
   type        = string
   description = "Prefix of S3 keys"
}

variable "asg_size" {
  type    = number
  default = 1
}

variable "security_group_ids" {
  type    = list(string)
  default = []
}

data "aws_ami" "amazon_linux_2" {
  most_recent = true

  filter {
    name   = "name"
    values = ["packer-amazon-linux-2-docker-*"]
  }

  owners = ["241748083751"] # Jeff Taylor
}

data "aws_availability_zones" "all" {}

data "aws_region" "current" {}

resource "aws_security_group" "worker" {
  name = "worker"
  dynamic "ingress" {
    for_each = var.allow_ssh ? [0] : []

    content {
      from_port   = 22
      to_port     = 22
      protocol    = "tcp"
      cidr_blocks = var.dev_cidrs
    }
  }

  # This rule is required for downloading docker images.
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
}

resource "aws_cloudwatch_log_group" "g" {
  name = var.image_name
}

module "iam_role" {
  source = "../iam_role"
  name = var.image_name
  policy_arns = ["arn:aws:iam::aws:policy/AmazonS3ReadOnlyAccess"]
}

resource "aws_iam_role_policy" "s3_write" {
  name = "S3WriteData-${var.image_name}"
  role = module.iam_role.role.id

  policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": "s3:*",
            "Resource": [
                "arn:aws:s3:::${var.s3_bucket}/${var.s3_prefix}/*"
            ]
        }
    ]
}
EOF
}

resource "aws_launch_configuration" "worker" {
  name_prefix     = "tf-${var.image_name}-"
  image_id        = data.aws_ami.amazon_linux_2.id
  instance_type   = var.instance_type
  security_groups = concat([aws_security_group.worker.id], var.security_group_ids)
  spot_price      = var.spot_price
  key_name        = "jeff@rosario"
  iam_instance_profile = module.iam_role.instance_profile.name

  root_block_device {
    volume_size = var.volume_size
  }

  # Launch configs can't be destroyed while attached to an ASG.
  lifecycle {
    create_before_destroy = true
  }

  user_data = <<-EOF
              #!/bin/bash
              docker run -d \
                --log-driver=awslogs \
                --log-opt awslogs-region="${data.aws_region.current.name}" \
                --log-opt awslogs-group="${aws_cloudwatch_log_group.g.name}" \
                --name ${var.image_name} \
                jefftaylor42/${var.image_name} \
                -u ${var.scheduler_dns}:${var.scheduler_port} \
                -k s3://${var.s3_bucket}/${var.s3_prefix}
              EOF
}

resource "aws_autoscaling_group" "worker" {
  # Forces a change in launch configuration to create a new ASG.
  name                 = "tf-asg-${aws_launch_configuration.worker.name}"

  launch_configuration = aws_launch_configuration.worker.id
  availability_zones   = data.aws_availability_zones.all.names

  min_size         = 0
  desired_capacity = var.up ? var.asg_size : 0
  max_size         = var.up ? var.asg_size : 0
}

