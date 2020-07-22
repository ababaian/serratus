#!/usr/bin/env python3

'''
Example:
aws s3 cp s3://serratus-rayan/sra_master_table.csv /tmp/
awk -F, 'NR>1 { print $1 }' /tmp/sra_master_table.csv | time python ./sra_metadata.py | tee /tmp/sra_disaster_table.tsv 


Output fields:
1. NCBI Nucleotide accession
2. NCBI Taxonomy  ID
3. Length of sequence
4. Family Taxon ID (if available)
5. Family Taxon ID (if available)

Note: sometimes a virus has unclear taxonomy, so there is no 'family' rank available in its lineage.

'''

from collections import defaultdict
import xml.etree.ElementTree as ET
import fileinput, time, sys
from Bio import Entrez

## Local Definitions:
Entrez.email = "serratus@me.tomeraltman.net"


## Read in all of the NCBI Nucelotide identifiers:
def ingest_ids():
    return [ line.strip().split(',')[0] for line in fileinput.input()]

def fetch_sra_xml_records():

    num_ids = []
    
    ## Use EPost to store a large number of IDs:
    id_list = ingest_ids()

    for i in range(0, divmod(len(id_list), 20)[0]+1):
        print(i, file=sys.stderr)
        id_sublist = id_list[(i*20):((i+1)*20-1)]
        if len(id_sublist) > 0:
            time.sleep(0.5)
            curr_num_ids = Entrez.read(Entrez.esearch(db="sra", term=" OR ".join(id_sublist)))['IdList']
            num_ids.extend(curr_num_ids)

        
    ## We need to convert the accessions to the internal numerical IDs used by EPost:
    search_results = Entrez.read(Entrez.epost("sra", id=",".join(num_ids)))
    query_key = search_results["QueryKey"]
    webenv = search_results["WebEnv"]
    records = []
    
    ## Download the GenBank data in XML:
    fetch_handle = Entrez.efetch(db="sra", webenv=webenv, query_key=query_key, retmode="xml")
    root = ET.fromstring(fetch_handle.read())
    
    
    return(root)


if __name__ == "__main__":
    hosts_dict = {}
    min_host_occur = 3
    
    xml_records = fetch_sra_xml_records()

    ## Print accessions:
    for exp_pkg in xml_records:
        try:
            acc = exp_pkg.find('.//RUN').attrib['accession']
            platform = exp_pkg.find('.//PLATFORM')[0].tag
            instrument_model = exp_pkg.find('.//PLATFORM')[0][0].text
            sample_taxon = exp_pkg.find('.//SAMPLE_NAME/TAXON_ID').text
            sample_sciname = exp_pkg.find('.//SAMPLE_NAME/SCIENTIFIC_NAME').text
            bioproject = exp_pkg.find('.//EXTERNAL_ID[@namespace="BioProject"]').text
            for run in exp_pkg.findall('.//RUN'):
                acc = run.attrib['accession']
                print(acc, bioproject, platform, instrument_model, sample_taxon, sample_sciname, sep="\t")
        except:
            print(acc, "Error processing SRA run", sep="\t")
