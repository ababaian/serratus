variable "name" { 
  type = string
}

variable "policy_arns" {
  type = set(string)
  default = []
}

resource "aws_iam_role" "role" {
  name = "SerratusIamRole-${var.name}"

  assume_role_policy = <<EOF
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

resource "aws_iam_role_policy" "cloudwatch" {
  name = "CloudWatchLogsCreate-${var.name}"
  role = aws_iam_role.role.id

  policy = <<EOF
{
	"Version": "2012-10-17",
	"Statement": [
		{
			"Action": [
				"logs:CreateLogStream",
				"logs:PutLogEvents"
			],
			"Effect": "Allow",
			"Resource": "*"
		}
	]
}
EOF
}

resource "aws_iam_instance_profile" "profile" {
  name = "profile-${var.name}"
  role = aws_iam_role.role.name
}

resource "aws_iam_role_policy_attachment" "attachment" {
  for_each = var.policy_arns

  role       = aws_iam_role.role.name
  policy_arn = each.value
}

output "role" {
  value = aws_iam_role.role
}

output "instance_profile" {
  value = aws_iam_instance_profile.profile
}
