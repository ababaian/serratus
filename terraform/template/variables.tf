///////////////////////////////////////////////////////////
// CLUSTER - Variables (TEMPLATE)
///////////////////////////////////////////////////////////

// AWS User Settings ======================================
// Account settings and security key
// (Note: also configure `main.tf` S3 backend to match)

variable "aws_account" {
  type = string
  default = "797308887321"
  description = "AWS Account ID"
}

variable "aws_profile" {
  type = string
  default = "default"
  description = "AWS profile [default] in most cases"
}

variable "aws_key_name" {
  type = string
  default   = "biof"
  description = "Public key name stored locally for SSH into instance"
}

variable "aws_public_key" {
  //to make: ssh-keygen -y -f /path_to_private_key.pem
  type = string
  default = "ssh-rsa <KEYGEN HERE>"
  description = "Public key on local system"
}

// AWS Pipeline Settings ==================================
// Project specific settings with which to initialize serratus

variable "aws_region" {
  type = string
  default = "us-east-1"
  description = "AWS region in which to run pipeline. NOTE: ami-id are region-specific."
}

variable "aws_bucket" {
  type = string
  default = "serratus-live"
  description = "S3 Bucket to use as base for pipeline"
}

variable "aws_subnet" {
  type = string
  default = "subnet-e6f6e8ab"
  description = "Externally defined AWS subnet."
}

variable "serratus-base-ami" {
  type = string
  default = "<AMI-ID>"
  description = "Packer generated ami to use as base for serratus"
}
