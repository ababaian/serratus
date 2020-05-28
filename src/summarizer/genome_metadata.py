#!/usr/bin/env python3

'''
Example:
time python ./genome_metadata.py < ../../flom2.fa.fai

Output fields:
1. NCBI Nucleotide accession
2. NCBI Taxonomy  ID
3. Length of sequence
4. Family Taxon ID (if available)
5. Family Taxon ID (if available)

Note: sometimes a virus has unclear taxonomy, so there is no 'family' rank available in its lineage.

'''

from collections import defaultdict
import fileinput
from Bio import Entrez

## Local Definitions:
Entrez.email = "serratus@me.tomeraltman.net"
seq2lengths = defaultdict(int)

## Read in all of the NCBI Nucelotide identifiers:
def ingest_ids():
    for line in fileinput.input():
        fields = line.strip().split('\t')
        seq2lengths[fields[0]] = fields[1]

def fetch_nucleotide_genbank_records():

    ## Use EPost to store a large number of IDs:
    id_list = list(seq2lengths.keys())
    search_results = Entrez.read(Entrez.epost("nucleotide", id=",".join(id_list)))
    query_key = search_results["QueryKey"]
    webenv = search_results["WebEnv"]

    ## Download the GenBank data in XML:
    fetch_handle = Entrez.efetch(db="nuccore", webenv=webenv, query_key=query_key, retmode="xml")

    ## Parse the XML into records:
    records = Entrez.read(fetch_handle)
    fetch_handle.close()
    
    return(records)

## We need to fetch:
## * TaxonID
## * Family name
## * TaxonID of family name:

def scrape_records(records):
    record_data = []
    for record in records:
        accession = record['GBSeq_accession-version']
        taxon_id = [ qual['GBQualifier_value'].split(':')[1]
                     for qual in record['GBSeq_feature-table'][0]['GBFeature_quals']
                     if qual['GBQualifier_name'] == 'db_xref' and qual['GBQualifier_value'].startswith('taxon:') ]
        record_data.append([accession, taxon_id[0], seq2lengths[accession]])
    return(record_data)

def fetch_taxon_family_data(record_data):
    taxon2family_data = {}
    
    taxon_ids = [ record[1] for record in record_data]
    search_results = Entrez.read(Entrez.epost("taxonomy", id=",".join(taxon_ids)))
    query_key = search_results["QueryKey"]
    webenv = search_results["WebEnv"]

    ## Download the Taxonomy data in XML:
    fetch_handle = Entrez.efetch(db="taxonomy", webenv=webenv, query_key=query_key, retmode="xml")
    tax_records = Entrez.read(fetch_handle)

    for taxon in tax_records:
        taxon_id = taxon['TaxId']
        family_data = [[lineage_taxon['TaxId'], lineage_taxon['ScientificName']] 
                       for lineage_taxon in taxon['LineageEx'] 
                       if lineage_taxon['Rank'] == 'family' ]
        if len(family_data) > 0:
            taxon2family_data[taxon_id] = family_data[0]
    
    return(taxon2family_data)


if __name__ == "__main__":

    ingest_ids()
    gb_records = fetch_nucleotide_genbank_records()
    record_data = scrape_records(gb_records)
    taxon2family_data = fetch_taxon_family_data(record_data)
    #print(taxon2family_data)
    for record in record_data:
        #print(record)
        if record[1] in taxon2family_data:
            print('\t'.join(record + taxon2family_data[record[1]]))
        else:
            print('\t'.join(record + ['', '']))

    exit()
