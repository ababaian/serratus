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

variable "instance_type" {
  type = string
  default = "t3.nano"
}

data "aws_ami" "ecs" {
  most_recent = true
  owners      = ["591542846629"] # Amazon

  filter {
    name   = "name"
    values = ["amzn2-ami-ecs-hvm-*"]
  }

  filter {
    name   = "architecture"
    values = ["x86_64"]
  }
}

// RESOURCES ==============================================

# Give our instance the set of permissions required to act as an ECS node.
resource "aws_iam_instance_profile" "monitor" {
  name = "profile-serratus-monitor"
  role = aws_iam_role.instance_role.name
}

resource aws_instance "monitor" {
  ami = data.aws_ami.ecs.id
  instance_type = var.instance_type
  instance_initiated_shutdown_behavior = "terminate"
  vpc_security_group_ids               = var.security_group_ids
  key_name                             = var.key_name
  iam_instance_profile                 = aws_iam_instance_profile.monitor.name

  tags = {
    "project": "serratus"
    "component": "serratus-monitor"
  }

  user_data = <<-EOF
    #cloud-config
    write_files:
      - path: /etc/ecs/ecs.config
        content: |
          ECS_CLUSTER=${aws_ecs_cluster.monitor.name}
          ECS_ENABLE_TASK_IAM_ROLE_NETWORK_HOST=true
    # Download and run node exporter, so this instance can monitor itself.
    # For the other, we'll use a custom AMI which has the exporter built
    # in.
    bootcmd:
      - [ "curl", "-sL", "https://github.com/prometheus/node_exporter/releases/download/v1.0.0-rc.0/node_exporter-1.0.0-rc.0.linux-amd64.tar.gz", "-o", "node_exporter.tgz" ]
      - [ "tar", "-xvzf", "node_exporter.tgz" ]
      - [ "systemd-run", "node_exporter-1.0.0-rc.0.linux-amd64/node_exporter" ]
  EOF
}

resource "aws_eip" "monitor" {
  instance = aws_instance.monitor.id
  vpc      = true

  tags = {
    "project": "serratus"
    "component": "serratus-monitor"
  }
}

resource aws_ecs_cluster "monitor" {
  name = "serratus-monitor"
}

resource "aws_iam_role" "instance_role" {
  name = "SerratusEcsInstanceRole"

  assume_role_policy = <<-EOF
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

# Allows this EC2 instance to act as an ECS agent
resource "aws_iam_role_policy_attachment" "instance_attachment" {
  role       = aws_iam_role.instance_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

resource "aws_iam_role" "task_role" {
  name = "SerratusIamRole-monitor"

  assume_role_policy = <<-EOF
    {
      "Version": "2012-10-17",
      "Statement": [
        {
          "Action": "sts:AssumeRole",
          "Principal": {
            "Service": "ecs-tasks.amazonaws.com"
          },
          "Effect": "Allow",
          "Sid": ""
        }
      ]
    }
  EOF
}

# Also allow the container to load cloudwatch metrics
resource "aws_iam_role_policy" "cloudwatch" {
  name = "CloudwatchGetMetrics"
  role = aws_iam_role.task_role.id

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

# Provides read-only access to EC2 "Describe" functions, ELB, Autoscaling,
# and Cloudwatch
resource "aws_iam_role_policy_attachment" "attachment" {
  role       = aws_iam_role.task_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonEC2ReadOnlyAccess"
}

resource aws_ecs_task_definition "monitor" {
  family = "monitor"
  container_definitions = templatefile("../monitoring/monitor-task-definition.json",
                                       { sched_ip = var.scheduler_ip})
  task_role_arn = aws_iam_role.task_role.arn
  network_mode = "host"

  volume {
    name = "prometheus-data"

    docker_volume_configuration {
      scope = "shared"
      autoprovision = true
      driver = "local"
    }
  }
}

resource aws_ecs_service "monitor" {
  name = "serratus-monitor"
  cluster = aws_ecs_cluster.monitor.id
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

