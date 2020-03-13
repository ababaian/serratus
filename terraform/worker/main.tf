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

#variable "scheduler_ip" {
#  description = "IP Address of the scheduler"
#  type        = string
#}

variable "instance_type" {
  description = "Type of node to use for the workers"
  type        = string
  default     = "t3.nano"
}

variable "spot_price" {
  type    = number
  default = 0.0016
}

variable "allow_ssh" {
  description = "Allow SSH access to the nodes"
  type        = bool
  default     = true
}

variable "my_ip" {
  description = "My IP Address, for SSH, HTTP, etc access"
  type        = string
}

variable "name" {
  description = "Default name of sub-components"
  type        = string
  default     = "worker"
}

variable "asg_size" {
  type    = number
  default = 1
}

data "aws_ami" "amazon_linux_2" {
  # A simple AMI, built from Amazon Linux 2, plus the following script:
  # yum update -yq
  # yum install docker -yq
  # systemctl enable docker
  #
  # TODO: Put this in packer, so we can switch Regions / Clouds more easily
  most_recent = true

  filter {
    name   = "name"
    values = ["amazon_linux_2_docker"]
  }

  owners = ["241748083751"] # Jeff Taylor
}

data "aws_availability_zones" "all" {}

resource "aws_security_group" "worker" {
  name = "worker"
  dynamic "ingress" {
    for_each = var.allow_ssh ? [0] : []

    content {
      from_port   = 22
      to_port     = 22
      protocol    = "tcp"
      cidr_blocks = ["${var.my_ip}/32"]
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

resource "aws_cloudwatch_log_group" "worker" {
  name = "worker"
}

resource "aws_launch_configuration" "worker" {
  image_id        = data.aws_ami.amazon_linux_2.id
  instance_type   = var.instance_type
  security_groups = [aws_security_group.worker.id]
  spot_price      = var.spot_price
}

resource "aws_autoscaling_group" "worker" {
  launch_configuration = aws_launch_configuration.worker.id
  availability_zones   = data.aws_availability_zones.all.names

  min_size         = 0
  desired_capacity = var.up ? var.asg_size : 0
  max_size         = var.up ? var.asg_size : 0
}

