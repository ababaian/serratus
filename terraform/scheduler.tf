provider "aws" {
  profile = "default"
  region  = "us-west-2"
}

variable "port" {
  description = "HTTP port to use for the scheduler"
  type        = number
  default     = 8000
}

variable "vpc" {
  description = "AWS VPC id to use"
  type        = string
}

data "aws_availability_zones" "all" {}
data "aws_subnet_ids" "all" {
  vpc_id = var.vpc_id
}

# Security Group
resource "aws_security_group" "scheduler" {
  # TODO: Allow HTTP traffic from job runners
  name = "serratus-scheduler"

  ingress {
    from_port = var.port
    to_port = var.port
    protocol = "tcp"
    # TODO: only from workers
    cidr_blocks = ["0.0.0.0/0"]
  }
}

# TODO
# * Create EC2 instance + AMI to fetch and run the container
# * Image file
# * Logging Driver
# * VPC, Subnet
# * IAM Role
# * Task Definition
# * Cluster

resource "aws_autoscaling_group" "example" {
  launch_configuration = aws_launch_configuration.example.id
  availability_zones   = data.aws_availability_zones.all.names

  min_size = 2
  max_size = 3

  load_balancers    = [aws_elb.example.name]
  health_check_type = "ELB"

  tag {
    key                 = "Name"
    value               = "terraform-asg-example"
    propagate_at_launch = true
  }
}

resource "aws_ecs_cluster" "scheduler" {
  name = "serratus-scheduler"
  capacity_providers = [aws_ecs_capacity_provider.scheduler]
}

resource "aws_ecs_capacity_provider" "scheduler" {
  name = "serratus-scheduler"
  auto_scaling_group_provider {
    auto_scaling_group_provider_arn = aws_autoscaling_group.scheduler.arn
  }
}

# ECS Service
resource "aws_ecs_service" "scheduler" {
  name = "serratus-scheduler"
  cluster = aws_ecs_cluster.serratus_cluster.id
  task_definition = aws_ecs_task_definition.serratus_task.arn
  desired_count = 1
  iam_role = aws_iam_role.serratus_scheduler_role.arn
  depends_on      = [aws_iam_role_policy.serratus_scheduler_role]
  launch_type = "EC2"

  capacity_provider_strategy {
    capacity_provider = aws_ecs_capacity_provider.serratus_scheduler_provider.arn
    weight = 1
  }

  network_configuration {
    subnets = data.aws_subnet_ids.all.id
    security_groups = [aws_security_group.scheduler.id]
    assign_public_ip = true # Required to pull the image, TODO see if we can turn this off.
  }
}
