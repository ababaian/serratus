///////////////////////////////////////////////////////////
// MAIN.TF - CLUSTER (TEMPLATE)
///////////////////////////////////////////////////////////
// Variables disallowed
terraform {
  required_version = "= 0.12.21"
  backend "s3" {
    profile = "default"
    region  = "us-east-1"
    bucket  = "<BUCKET-NAME>"
    key     = "<PATH_TO_/key.pem>"
  }
}

// PROVIDER ===============================================
provider "aws" {
  version             = "~> 2.49"
  region              = var.aws_region
  profile             = var.aws_profile
  allowed_account_ids = ["${var.aws_account}"]
}

// AMI ====================================================

resource "aws_instance" "base_batch" {
  key_name          = var.aws_key_name
  security_group    = "default"
  // Amazon Linux base : ami-079f731edfe27c29c
  source_ami_id     = "${var.serratus-base-ami}"
  source_ami_region = "${var.aws_region}"
  encrypted         = "true"
  instance_type     = "t2.micro"
  tags = {
    Name = "serrautus-base"
  }
}