provider "aws" {
  profile = "default"
  region  = "us-west-2"
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

data "aws_region" "current" {}

variable "scheduler_port" {
  description = "HTTP port to use for the scheduler"
  type        = number
  default     = 8000
}

variable "input_security_group_ids" {
  description = "Security groups to allow data in"
  type        = list(number)
  default     = []
}

variable "allow_ssh" {
  description = "Allow SSH access to the nodes"
  type        = bool
  default     = true
}

variable "up" {
  description = "Spin up instances"
  type = bool
  default = true
}

resource "aws_security_group" "scheduler" {
  name = "serratus-scheduler"
  ingress {
    from_port       = var.scheduler_port
    to_port         = var.scheduler_port
    security_groups = var.input_security_group_ids
    protocol        = "tcp"
    cidr_blocks     = ["0.0.0.0/0"] # TODO Testing only, delete me
  }
  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }
  dynamic "ingress" {
    for_each = var.allow_ssh ? [0] : []

    content {
      from_port   = 22
      to_port     = 22
      protocol    = "tcp"
      cidr_blocks = ["0.0.0.0/0"]
    }
  }
}

resource "aws_iam_role" "scheduler" {
  name = "scheduler"

  assume_role_policy = <<EOF
{
		"Version": "2012-10-17",
		"Statement": [
				{
						"Action": "sts:AssumeRole",
						"Principal": {
							"Service": "ec2.amazonaws.com"
						},
						"Effect": "Allow",
						"Sid": ""
				}
		]
}
EOF
}

resource "aws_iam_instance_profile" "scheduler" {
  name = "scheduler"
  role = aws_iam_role.scheduler.name
}

resource "aws_iam_role_policy" "scheduler" {
  name = "scheduler"
  role = aws_iam_role.scheduler.id

  policy = <<EOF
{
	"Version": "2012-10-17",
	"Statement": [
		{
			"Action": [
				"logs:CreateLogStream",
				"logs:PutLogEvents"
			],
			"Effect": "Allow",
			"Resource": "*"
		}
	]
}
EOF
}

resource "aws_cloudwatch_log_group" "scheduler" {
  name = "scheduler"
}

resource "aws_instance" "scheduler" {
  ami                                  = data.aws_ami.amazon_linux_2.id
  instance_initiated_shutdown_behavior = "terminate"
  instance_type                        = "t3.nano"
  vpc_security_group_ids               = [aws_security_group.scheduler.id]
  key_name                             = "jeff@rosario"
  iam_instance_profile                 = aws_iam_instance_profile.scheduler.name

  count = var.up ? 1 : 0

  user_data = <<-EOF
              #!/bin/bash
              docker run -d -p "${var.scheduler_port}":8000 \
                --log-driver=awslogs \
                --log-opt awslogs-region="${data.aws_region.current.name}" \
                --log-opt awslogs-group="${aws_cloudwatch_log_group.scheduler.name}" \
                --name sch \
                jefftaylor42/serratus-scheduler
              EOF

  credit_specification {
    cpu_credits = "standard"
  }
}

#output "scheduler_public_dns" {
#  value = aws_instance.scheduler.public_dns
#}
