test:


$(DATA_DIR)/cov-hits.gb:
	mkdir -p $(DATA_DIR)
	cd $(DATA_DIR)
	aws s3 cp s3://serratus-public/out/200505_zoonotic/summary/ERR2756788.summary .
	awk -F';' '/^acc=/ {split($$1,parts,"="); if(parts[2] != "pan_genome") print parts[2]}' \
		ERR2756788.summary \
	| epost -db nuccore \
	| efetch -format gb \
		> $@

$(DATA_DIR)/cov-hits-taxa.txt:
	mkdir -p $(DATA_DIR)
	cd $(DATA_DIR)
	aws s3 cp s3://serratus-public/out/200505_zoonotic/summary/ERR2756788.summary .
	awk -F';' '/^acc=/ {split($$1,parts,"="); if(parts[2] != "pan_genome") print parts[2]}' \
		ERR2756788.summary \
	| epost -db nuccore \
	| elink -target taxonomy \
	| efetch -mode xml
