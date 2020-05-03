///////////////////////////////////////////////////////////
// Worker Node
///////////////////////////////////////////////////////////
//
//

// VARIABLES ==============================================

variable "instance_type" {
  description = "Type of node to use for the workers"
  type        = string
}

variable "spot_price" {
  type = number
}

variable "volume_size" {
  type = number
  description = "Size of the root EBS volume in GB"
  default = 8
}

variable "volume_type" {
  type        = string
  description = "Type of the root EBS volume: io1, gp2, st1, sc"
  default     = "gp2"
}

variable "allow_ssh" {
  description = "Allow SSH access to the nodes"
  type        = bool
  default     = true
}

variable "scheduler" {
  type = string
}

variable "workers" {
  type        = number
  description = "Number of worker threads to use"
  default     = 1
}

variable "dev_cidrs" {
  description = "Remote IP Address, for SSH, HTTP, etc access"
  type        = set(string)
}

variable "key_name" {
  description = "Name of the AWS key pair to assign instances"
  type        = string
}

variable "dockerhub_account" {
  description = "Docker Hub account to pull from"
  type        = string
}

variable "image_name" {
  description = "Docker image to run once the container is started"
  type        = string
}

variable "s3_bucket" {
  type        = string
  description = "Name of the S3 bucket to store temporary data"
}

variable "s3_prefix" {
  type        = string
  description = "Prefix of S3 keys"
}

variable "s3_delete_prefix" {
  type        = string
  description = "Prefix where instance should be able to delete objects"
  default     = ""
}

variable "desired_size" {
  description = "Desired size of the ASG"
  type        = number
  default     = 0
}

variable "max_size" {
  description = "Desired size of the ASG"
  type        = number
  default     = 1
}

variable "security_group_ids" {
  type    = list(string)
  default = []
}

variable "options" {
  type = string
  description = "Options to pass to the container"
}

data "aws_ami" "amazon_linux_2" {
  most_recent = true
  owners      = ["self"]

  filter {
    name   = "name"
    values = ["packer-amazon-linux-2-docker-*"]
  }
}

data "aws_availability_zones" "all" {}

data "aws_region" "current" {}

// RESOURCES ==============================================

resource "aws_cloudwatch_log_group" "g" {
  name = var.image_name
}

module "iam_role" {
  source      = "../iam_role"
  name        = var.image_name
  policy_arns = ["arn:aws:iam::aws:policy/AmazonS3ReadOnlyAccess"]
}


resource "aws_iam_role_policy" "ec2Describe" {
  name = "DescribeEC2Instances-${var.image_name}"
  role = module.iam_role.role.id

  policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Action": [
        "ec2:Describe*"
      ],
      "Effect": "Allow",
      "Resource": "*"
    }
  ]
}
EOF
}

resource "aws_iam_role_policy" "ec2Terminate" {
  name = "TerminateEC2Instances-${var.image_name}"
  role = module.iam_role.role.id

  policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Action": [
        "ec2:Terminate*"
      ],
      "Effect": "Allow",
      "Resource": "*"
    }
  ]
}
EOF
}

resource "aws_iam_role_policy" "AdjustAutoScaling" {
  name = "AdjustAutoScaling-${var.image_name}"
  role = module.iam_role.role.id

  policy = <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Action": [
        "autoscaling:*"
      ],
      "Effect": "Allow",
      "Resource": "*"
    }
  ]
}
EOF
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

resource "aws_iam_role_policy" "s3_delete" {
  name = "S3DeleteData-${var.image_name}"
  role = module.iam_role.role.id
  count = var.s3_delete_prefix != "" ? 1 : 0

  policy = <<EOF
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": "s3:DeleteObject",
            "Resource": [
                "arn:aws:s3:::${var.s3_bucket}/${var.s3_delete_prefix}/*"
            ]
        }
    ]
}
EOF
}

resource "aws_launch_configuration" "worker" {
  name_prefix          = "${var.image_name}-"
  image_id             = data.aws_ami.amazon_linux_2.id
  instance_type        = var.instance_type
  security_groups      = var.security_group_ids
  spot_price           = var.spot_price
  key_name             = var.key_name
  iam_instance_profile = module.iam_role.instance_profile.name
  enable_monitoring    = false

  root_block_device {
    volume_size = var.volume_size
    volume_type = var.volume_type
  }

  # Launch configs can't be destroyed while attached to an ASG.
  lifecycle {
    create_before_destroy = true
  }

  user_data = <<-EOF
              #!/bin/bash
              export 
              instance_id=$(curl -s http://169.254.169.254/latest/meta-data/instance-id)
              hostname ${var.image_name}-$instance_id
              docker run -d \
                --log-driver=awslogs \
                --log-opt awslogs-region="${data.aws_region.current.name}" \
                --log-opt awslogs-group="${aws_cloudwatch_log_group.g.name}" \
                --log-opt awslogs-stream="$instance_id" \
                --name ${var.image_name} \
                -e SCHEDULER=${var.scheduler} \
                ${var.dockerhub_account}/${var.image_name} \
                ${var.options}
              EOF
}


// TODO: COOLDOWN POLICY NOT ATTACHED TO GROUP
resource "aws_autoscaling_policy" "worker" {
  name = aws_launch_configuration.worker.name
  autoscaling_group_name = aws_autoscaling_group.worker.name
  scaling_adjustment     = 5
  adjustment_type        = "ChangeInCapacity"
  cooldown               = 30
}

resource "aws_autoscaling_group" "worker" {
  name = aws_launch_configuration.worker.name

  launch_configuration = aws_launch_configuration.worker.id
  availability_zones   = data.aws_availability_zones.all.names

  min_size         = 0
  desired_capacity = var.desired_size
  max_size         = var.max_size

  lifecycle {
    ignore_changes = [
      # This can be changed by scripts or by the scheduler (in the future)
      desired_capacity,
    ]
  }

  tag {
    key                 = "project"
    value               = "serratus"
    propagate_at_launch = true
  }

  tag {
    key                 = "component"
    value               = var.image_name
    propagate_at_launch = true
  }

  tag {
    key                 = "Name"
    value               = "${var.image_name}-instance"
    propagate_at_launch = true
  }
}

output "asg_name" {
  value = aws_autoscaling_group.worker.name
}

