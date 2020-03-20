resource "aws_s3_bucket" "terraform_state" {
  bucket = "serratus-tf-state"
  region = "us-west-2"

  # Enable versioning so we can see history of the state files
  versioning {
    enabled = true
  }

  grant {
    # jefft.second
    type = "CanonicalUser"
    id = "dcd86c2c4a99bfb6458a936948b17280a6120ebfa47491c638ae772477fefea0"
    permissions = ["FULL_CONTROL"]
  }

  server_side_encryption_configuration {
    rule {
      apply_server_side_encryption_by_default {
        sse_algorithm = "AES256"
      }
    }
  }
}

resource "aws_dynamodb_table" "terraform_locks" {
  name         = "serratus-tf-locks"
  billing_mode = "PAY_PER_REQUEST"
  hash_key     = "LockID"

  attribute {
    name = "LockID"
    type = "S"
  }
}
