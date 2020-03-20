variable "scheduler_port" {
  description = "HTTP port to use for the scheduler"
  type        = number
  default     = 8000
}

variable "instance_type" {
  type = string
}

variable "allow_ssh" {
  description = "Allow SSH access to the nodes"
  type        = bool
  default     = true
}

variable "dev_cidrs" {
  description = "Remote IP Address, for testing access"
  type        = set(string)
}

variable "security_groups" {
  type = list(string)
  default = []
}

variable "security_group_ids" {
  type = list(string)
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

data "aws_region" "current" {}

resource "aws_security_group" "scheduler" {
  name = "serratus-scheduler"
  ingress {
    from_port       = var.scheduler_port
    to_port         = var.scheduler_port
    protocol        = "tcp"
    cidr_blocks     = var.dev_cidrs
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
      cidr_blocks = var.dev_cidrs
    }
  }
}

module "iam_role" {
  source = "../iam_role"
  name = "scheduler"
}

resource "aws_cloudwatch_log_group" "scheduler" {
  name = "scheduler"
}

resource "aws_eip" "sch" {
  instance = aws_instance.scheduler.id
  vpc = true
}

resource "aws_instance" "scheduler" {
  ami                                  = data.aws_ami.amazon_linux_2.id
  instance_initiated_shutdown_behavior = "terminate"
  instance_type                        = var.instance_type
  vpc_security_group_ids               = concat([aws_security_group.scheduler.id], var.security_group_ids)
  key_name                             = "jeff@rosario"
  iam_instance_profile                 = module.iam_role.instance_profile.name
  monitoring                           = true

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

output "public_dns" {
  value = aws_eip.sch.public_dns
}

