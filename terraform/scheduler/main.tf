///////////////////////////////////////////////////////////
// SCHEDULER MODULE
///////////////////////////////////////////////////////////

provider random {}

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


// RESOURCES ##############################

module "ecs_cluster" {
  source = "../ecs_cluster"

  name = "scheduler"
  instance_type = var.instance_type
  security_group_ids = var.security_group_ids
  key_name = var.key_name
}

# Give the scheduler an Elastic-IP so we can destroy and recreate it without
# Also destroying everything that depends on its IP.
resource "aws_eip" "sch" {
  instance = module.ecs_cluster.instance.id
  vpc      = true

  tags = {
    "project": "serratus"
    "component": "serratus-scheduler"
  }
}

resource "aws_iam_role_policy" "scheduler" {
  name = "DescribeInstances-scheduler"
  role = module.ecs_cluster.task_role.id

  policy = <<EOF
{
	"Version": "2012-10-17",
	"Statement": [
		{
			"Action": [
        "ec2:DescribeInstances",
        "autoscaling:SetDesiredCapacity",
        "autoscaling:DescribeAutoScalingGroups"
			],
			"Effect": "Allow",
			"Resource": "*"
		}
	]
}
EOF
}

resource "aws_cloudwatch_log_group" "scheduler" {
  name = "serratus-scheduler"
}


resource "random_password" "pg_password" {
  length = 16
  special = false
}

resource aws_ecs_task_definition "scheduler" {
  family = "scheduler"
  container_definitions = templatefile("../scheduler/scheduler-task-definition.json", {
    dockerhub_account  = var.dockerhub_account
    sched_port = var.scheduler_port,
    aws_region = data.aws_region.current.name
    log_group  = aws_cloudwatch_log_group.scheduler.name
    pg_password = random_password.pg_password.result
  })
  task_role_arn = module.ecs_cluster.task_role.arn
  network_mode = "host"

  volume {
    name = "postgres-data"

    docker_volume_configuration {
      scope = "shared"
      autoprovision = true
      driver = "local"
    }
  }
}

resource aws_ecs_service "scheduler" {
  name = "serratus-scheduler"
  cluster = module.ecs_cluster.cluster.id
  task_definition = aws_ecs_task_definition.scheduler.arn

  # TODO: Allow this to scale-out to many instances.
  scheduling_strategy = "DAEMON"
}

// OUTPUT ##############################

output "private_ip" {
  value = module.ecs_cluster.instance.private_ip
}

output "public_ip" {
  value = aws_eip.sch.public_ip
}

output "public_dns" {
  value = aws_eip.sch.public_dns
}

output "pg_password" {
  value = random_password.pg_password.result
}

