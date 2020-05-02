///////////////////////////////////////////////////////////
// SCHEDULER MODULE
///////////////////////////////////////////////////////////
// VARIABLES ##############################

variable "scheduler_port" {
  description = "HTTP port to use for the scheduler"
  type        = number
  default     = 8000
}

variable "instance_type" {
  type = string
}

variable "key_name" {
  description = "Name of the AWS key pair to assign instances"
  type        = string
}

variable "dockerhub_account" {
  description = "Docker Hub account to pull from"
  type        = string
}

variable "security_group_ids" {
  type    = list(string)
  default = []
}

data "aws_ami" "amazon_linux_2" {
  most_recent = true
  owners      = ["self"]

  filter {
    name   = "name"
    values = ["packer-amazon-linux-2-docker-*"]
  }
}

data "aws_region" "current" {}

module "iam_role" {
  source = "../iam_role"
  name   = "scheduler"
}

// RESOURCES ##############################

resource "aws_cloudwatch_log_group" "scheduler" {
  name = "scheduler"
}

# Give the scheduler an Elastic-IP so we can destroy and recreate it without
# Also destroying everything that depends on its IP.
resource "aws_eip" "sch" {
  instance = aws_instance.scheduler.id
  vpc      = true

  tags = {
    "project": "serratus"
    "component": "serratus-scheduler"
  }
}

resource "aws_instance" "scheduler" {
  ami                                  = data.aws_ami.amazon_linux_2.id
  instance_initiated_shutdown_behavior = "terminate"
  instance_type                        = var.instance_type
  vpc_security_group_ids               = var.security_group_ids
  key_name                             = var.key_name
  iam_instance_profile                 = module.iam_role.instance_profile.name
  monitoring                           = false

  user_data = <<-EOF
              #!/bin/bash
              instance_id=$(curl -s http://169.254.169.254/latest/meta-data/instance-id)
              hostname serratus-scheduler
              docker run -d -p "${var.scheduler_port}":8000 \
                --log-driver=awslogs \
                --log-opt awslogs-region="${data.aws_region.current.name}" \
                --log-opt awslogs-group="${aws_cloudwatch_log_group.scheduler.name}" \
                --log-opt awslogs-stream="$instance_id" \
                --name sch \
                ${var.dockerhub_account}/serratus-scheduler
              EOF

  tags = {
    "Name": "serratus-scheduler"
    "project": "serratus"
    "component": "serratus-scheduler"
  }
}

// OUTPUT ##############################

output "private_ip" {
  value = aws_instance.scheduler.private_ip
}

output "public_ip" {
  value = aws_eip.sch.public_ip
}

output "public_dns" {
  value = aws_eip.sch.public_dns
}

