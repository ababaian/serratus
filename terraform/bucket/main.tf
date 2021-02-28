///////////////////////////////////////////////////////////
// S3 WORK BUCKET MODULE
///////////////////////////////////////////////////////////
// VARIABLES ==============================================

variable "prefixes" {
  type        = set(string)
  default     = []
  description = "Prefixes to monitor"
}

resource "aws_s3_bucket" "work" {
  bucket_prefix = "tf-serratus-work-"
  force_destroy = true

  tags = {
    "project" : "serratus"
    "component" : "serratus-scheduler"
  }
}

// RESOURCES ==============================================

resource "aws_s3_bucket_metric" "full" {
  bucket = aws_s3_bucket.work.bucket
  name   = "full"
}

resource "aws_s3_bucket_metric" "prefix" {
  for_each = var.prefixes
  bucket   = aws_s3_bucket.work.bucket
  name     = "prefix-${each.value}"

  filter {
    prefix = each.value
  }
}

// OUTPUT =================================================

output "name" {
  value = aws_s3_bucket.work.bucket
}
