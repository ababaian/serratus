

#### Assumptions:
## * sudo apt-get install make
## * AWS CLI configured

#### Make configuration:

## Use Bash as default shell, and in strict mode:
SHELL := /bin/bash
.SHELLFLAGS = -ec

## If the parent env doesn't ste TMPDIR, do it ourselves:
TMPDIR ?= /tmp

## This makes all recipe lines execute within a shared shell process:
## https://www.gnu.org/software/make/manual/html_node/One-Shell.html#One-Shell
.ONESHELL:

## If a recipe contains an error, delete the target:
## https://www.gnu.org/software/make/manual/html_node/Special-Targets.html#Special-Targets
.DELETE_ON_ERROR:

## This is necessary to make sure that these intermediate files aren't clobbered:
.SECONDARY:


#### Targets:

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

$(DATA_DIR)/sumzer.tsv:
	mkdir -p $(DATA_DIR)
	cd $(DATA_DIR)
	aws s3 cp s3://serratus-public/seq/flom2/flom2.fa.fai .
	awk 'BEGIN{OFS="\t"}{ print $$1, $$2 }' \
		$(DATA_DIR)/flom2.fa.fai \
		> acc-len.txt
	epost -db nuccore < acc.tmp  \
		| efetch -mode xml \
		| tee nuccore.xml \
		| xtract -pattern GBSeq \
			-tab "\n" -element GBSeq_accession-version \
			-block GBQualifier -tab "\n" \
			-element GBQualifier_name,GBQualifier_value \
		| awk 'BEGIN{OFS="\t"}\
			$$2=="" && check==0 { acc=$$1; check=1; next }\
			$$2~/^taxon/ { gsub(/^taxon:/,"",$$2); check=0; print acc, $$2 }' \
		| sort | uniq \
		> acc-taxonid.txt
	epost -db nuccore < acc.tmp \
		| elink -db nuccore -target taxonomy \
		| efetch -mode xml \
		| tee flom2.xml \
		| sed 's/<TaxaSet><Taxon>/<TaxaSet>\n<Taxon>/' \
		| awk 'BEGIN{FS="[<>]";OFS="\t"} \
			/^<Taxon/ { start=1; next } \
			/TaxId/ && start == 1 { taxonid = $$3; start=0; next } \
			/TaxId/ && start == 0 { rank_taxonid = $$3; next } \
			/Scientific/ && start == 0 { putative_name=$$3; next} \
			$$3=="family" {print taxonid, putative_name, rank_taxonid; next}' \
		> taxonid-family.txt
	awk 'BEGIN{OFS="\t" } \
		NR==FNR { len[$$1] = $$2; next } \
		{ print $$1, $$2, len[$$1] }'\
		acc-len.txt \
		acc-taxonid.txt \
		| awk 'BEGIN{OFS="\t"} \
		NR == FNR { taxon_family[$$1] = $$2; \
                            taxon_rank_taxon[$$1] = $$3; \
                            next }\
		{print $$1, $$2, $$3, \
			taxon_family[$$2], \
			taxon_rank_taxon[$$2]}' \
		taxonid-family.txt \
		- \
		> sumzer.tsv


