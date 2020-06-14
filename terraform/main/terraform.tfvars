///////////////////////////////////////////////////////////
// Serratus - Global Variables
///////////////////////////////////////////////////////////

// AWS User Settings ##############################
key_name          = "serratus"

// ensure your IAM has permissions
// to write to the output bucket
// (no terminal slash)
output_bucket     = "s3://ur_bucket_name/out/yymmdd_ii"

// Local Settings ##############################
// <Your IP>/32. Use `curl ipecho.net/plain; echo`
dev_cidrs         = ["173.181.15.58/32", "75.155.242.64/32"]

// Docker Hub ##############################
// Account to upload container images to or use defaults ["serratusbio"]
dockerhub_account = "serratusbio"