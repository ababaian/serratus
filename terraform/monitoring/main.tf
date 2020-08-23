///////////////////////////////////////////////////////////
// MONITORING MODULE
///////////////////////////////////////////////////////////

// VARIABLES ==============================================
variable "key_name" {
  description = "Name of the AWS key pair to assign instances"
  type        = string
}

variable "security_group_ids" {
  type    = list(string)
  default = ["sg-05978c64cc44e9807", "sg-de40aaf6"] # rosario, default
}

variable "scheduler_ip" {
  type = string
}

variable "metrics_ip" {
  type = string
}

variable "instance_type" {
  type    = string
  default = "t3.nano"
}

variable "dockerhub_account" {
  description = "Docker Hub account to pull from"
  type        = string
}

data "aws_region" "current" {}

// RESOURCES ==============================================

module "ecs_cluster" {
  source = "../ecs_cluster"

  name               = "monitor"
  instance_type      = var.instance_type
  security_group_ids = var.security_group_ids
  key_name           = var.key_name
}

resource "aws_eip" "monitor" {
  instance = module.ecs_cluster.instance.id
  vpc      = true

  tags = {
    "project" : "serratus"
    "component" : "serratus-scheduler"
  }
}

# Provides read-only access to EC2 "Describe" functions, ELB, Autoscaling,
# and Cloudwatch
resource "aws_iam_role_policy_attachment" "attachment" {
  role       = module.ecs_cluster.task_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonEC2ReadOnlyAccess"
}

# Also allow the container to load cloudwatch metrics
resource "aws_iam_role_policy" "cloudwatch" {
  name = "CloudwatchGetMetrics"
  role = module.ecs_cluster.task_role.id

  policy = <<-EOF
    {
      "Version": "2012-10-17",
      "Statement": [
        {
          "Action": [
            "cloudwatch:GetMetricData"
          ],
          "Effect": "Allow",
          "Resource": "*"
        }
      ]
    }
  EOF
}

resource "aws_cloudwatch_log_group" "g" {
  name = "serratus-monitor"
}


resource aws_ecs_task_definition "monitor" {
  family = "monitor"
  container_definitions = templatefile("../monitoring/monitor-task-definition.json", {
    dockerhub_account = var.dockerhub_account
    sched_ip          = var.scheduler_ip,
    aws_region        = data.aws_region.current.name
    metrics_ip        = var.metrics_ip,
  })
  task_role_arn = module.ecs_cluster.task_role.arn
  network_mode  = "host"

  volume {
    name = "prometheus-data"

    docker_volume_configuration {
      scope         = "shared"
      autoprovision = true
      driver        = "local"
    }
  }
}

resource aws_ecs_service "monitor" {
  name            = "serratus-monitor"
  cluster         = module.ecs_cluster.cluster.id
  task_definition = aws_ecs_task_definition.monitor.arn

  # "Daemon" indicates that only one can be running at a time.  If we use the
  # default "replica" mode, ECS tries to create a new container before
  # destroying the old one, which fails because the ports are already taken.
  scheduling_strategy = "DAEMON"
}

// OUTPUT =================================================

output "private_ip" {
  value = aws_eip.monitor.private_ip
}

output "public_dns" {
  value = aws_eip.monitor.public_dns
}

