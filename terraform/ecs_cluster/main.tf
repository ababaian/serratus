//----------------------------------
// Variables
//----------------------------------

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

variable name {
  type = string
}

variable instance_type {
  type = string
}

variable security_group_ids {
  type    = list(string)
  default = []
}

variable "key_name" {
  description = "Name of the AWS key pair to assign instances"
  type        = string
}

//----------------------------------
// Resources
//----------------------------------

resource "aws_iam_role" "instance_role" {
  name = "SerratusEcsInstanceRole-${var.name}"

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

# Give our instance the set of permissions required to act as an ECS node.
resource "aws_iam_instance_profile" "p" {
  name = "instance-profile-serratus-${var.name}"
  role = aws_iam_role.instance_role.name
}

resource "aws_iam_role" "task_role" {
  name = "SerratusIamRole-${var.name}"

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


# Allows this EC2 instance to act as an ECS agent
resource "aws_iam_role_policy_attachment" "instance_attachment" {
  role       = aws_iam_role.instance_role.name
  policy_arn = "arn:aws:iam::aws:policy/service-role/AmazonEC2ContainerServiceforEC2Role"
}

resource aws_instance "i" {
  ami = data.aws_ami.ecs.id
  instance_type = var.instance_type
  instance_initiated_shutdown_behavior = "terminate"
  vpc_security_group_ids               = var.security_group_ids
  key_name                             = var.key_name
  iam_instance_profile                 = aws_iam_instance_profile.p.name

  tags = {
    "Name": "serratus-${var.name}"
    "project": "serratus"
    "component": "serratus-${var.name}"
  }

  user_data = <<-EOF
    #cloud-config
    write_files:
      - path: /etc/ecs/ecs.config
        content: |
          ECS_CLUSTER=${aws_ecs_cluster.c.name}
          ECS_ENABLE_TASK_IAM_ROLE_NETWORK_HOST=true
    # Use Node Exporter for CPU / Network / etc metrics.
    bootcmd:
      - [ "curl", "-sL", "https://github.com/prometheus/node_exporter/releases/download/v1.0.0-rc.0/node_exporter-1.0.0-rc.0.linux-amd64.tar.gz", "-o", "node_exporter.tgz" ]
      - [ "tar", "-xvzf", "node_exporter.tgz" ]
      - [ "systemd-run", "node_exporter-1.0.0-rc.0.linux-amd64/node_exporter" ]
  EOF
}

resource aws_ecs_cluster "c" {
  name = "serratus-${var.name}"
}

//----------------------------------
// Outputs
//----------------------------------

output "instance" {
  value = aws_instance.i
}

output "task_role" {
  value = aws_iam_role.task_role
}

output "cluster" {
  value = aws_ecs_cluster.c
}
