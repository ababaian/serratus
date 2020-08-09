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
import fileinput, time
from Bio import Entrez

## Local Definitions:
Entrez.email = "serratus@me.tomeraltman.net"


## Read in all of the NCBI Nucelotide identifiers:
def ingest_ids():
    return [ line.strip().split('\t')[0] for line in fileinput.input()]

def fetch_nucleotide_genbank_records():

    ## Use EPost to store a large number of IDs:
    id_list = ingest_ids()
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

def get_host(record):
    for qual in record['GBSeq_feature-table'][0]['GBFeature_quals']:
        if qual['GBQualifier_name'] == 'host':
            return qual['GBQualifier_value']
    return None

def get_taxon(record):
    for qual in record['GBSeq_feature-table'][0]['GBFeature_quals']:
        if qual['GBQualifier_name'] == 'db_xref' and qual['GBQualifier_value'].startswith('taxon:'):
            return qual['GBQualifier_value'].split(':')[1]
    return None

## Get NCBI Taxonomy ID for host term, plus the NCBI Taxonomy ID for the subsuming Order:
def get_host_taxa(host_name):
    search_result = Entrez.read(Entrez.esearch(db="taxonomy", term=host_name, retmode="xml"))
    if int(search_result['Count']) == 0:
        return []
    if int(search_result['Count']) > 1:
        return "more than one!"
    tax_id = search_result['IdList'][0]
    tax_metadata = Entrez.read(Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml"))
    sci_name = tax_metadata[0]['ScientificName']
    rank = tax_metadata[0]['Rank']
    if rank == "order":
        order_id = tax_id
        order_name = sci_name
    else:
        for taxon_dict in tax_metadata[0]['LineageEx']:
            if taxon_dict['Rank'] == 'order':
                order_id = taxon_dict['TaxId']
                order_name = taxon_dict['ScientificName']
    return [tax_id, sci_name, order_id, order_name]
    


if __name__ == "__main__":
    hosts_dict = {}
    min_host_occur = 3
    
    gb_records = fetch_nucleotide_genbank_records()
    ## Find unique hosts:
    for record in gb_records:
        host_name = get_host(record)
        hosts_dict[host_name] = {}
        hosts_dict[host_name]['count'] = 0
    
    # exit()
    for record in gb_records:
        host_name = get_host(record)
        
        if host_name == "camel":        
            hosts_dict[host_name]['alt_name'] = "camelus"
        if host_name == "Mustela putorius (ferret)":
            hosts_dict[host_name]['alt_name'] = "ferret"
        if host_name == "Pipistrellus cf. hesperidus; specimen voucher: OTBA03-20130220":
            hosts_dict[host_name]['alt_name'] = "Pipistrellus hesperidus"
        if host_name in ["bats", "bat", "microbat", "bat BF_258I", "bat BF_506I"]:
            hosts_dict[host_name]['alt_name'] = "chiroptera"
        if host_name == "Rhinolophus ferrumequinum (horseshoe bat)":
            hosts_dict[host_name]['alt_name'] = "Rhinolophus ferrumequinum"
        if host_name == "sparrow":
            hosts_dict[host_name]['alt_name'] = "passeridae"
        if host_name == "white-eye":
            hosts_dict[host_name]['alt_name'] = "Zosterops japonicus"
        if host_name == "wigeon":
            hosts_dict[host_name]['alt_name'] = "widgeon"
        if host_name == "night-heron":
            hosts_dict[host_name]['alt_name'] = "black-crowned night-heron"
        if host_name == "magpie-robin":
            hosts_dict[host_name]['alt_name'] = "oriental magpie robin"
        if host_name in ["Homo sapiens; patient #2 with severe acute respiratory syndrome (SARS)", "Homo sapiens; sex:M; age:34Y"]:
            hosts_dict[host_name]['alt_name'] = "Homo sapiens"
        if host_name == "rat":
            hosts_dict[host_name]['alt_name'] = "rattus"
        if host_name in ["piglet", "porcine", "Sus scrofa domesticus L.", "Porcine", "swine; suckling piglet", "newborn piglet", "swine; piglet"]:
            hosts_dict[host_name]['alt_name'] = "pig"
        if host_name == "bottlenose dolphin":
            hosts_dict[host_name]['alt_name'] = "Tursiops truncatus"
        if host_name in ["Bos taurus; breed: Holstein", "bovine; calf", "pre-weaned Korean native calf", "calf"]:
            hosts_dict[host_name]['alt_name'] = "Bos taurus"
        if host_name == "Canis lupus famaliaris":
            hosts_dict[host_name]['alt_name'] = "Canis lupus familiaris"
        if host_name == "feline":
            hosts_dict[host_name]['alt_name'] = "Felidae"
        if host_name == "Rhinolophus sp.":
            hosts_dict[host_name]['alt_name'] = "Rhinolophus"
        
        if host_name != None:
            ##print(host_name + " " + str(hosts_dict[host_name]['count']))
            hosts_dict[host_name]['count'] += 1        
        
    ## Iterate over hosts, get metadata:
    print('Num hosts: ' + str(len(list(hosts_dict.keys()))))
    for host in hosts_dict.keys():
        if hosts_dict[host]['count'] >= min_host_occur or 'alt_name' in hosts_dict[host]:
            time.sleep(1)
            if 'alt_name' in hosts_dict[host]:
                curr_host_taxa = get_host_taxa(hosts_dict[host]['alt_name'])
            else:
                curr_host_taxa = get_host_taxa(host)
            
            if curr_host_taxa == "more than one!":
                print("More than one: " + host)
            elif curr_host_taxa == []:
                print("No host found: " + host)
            else: 
                hosts_dict[host]['metadata'] = curr_host_taxa

                
    ## Print out final list:
    for record in gb_records:
        acc = record['GBSeq_locus']
        host = get_host(record)
        taxon = get_taxon(record)
        if taxon == None:
            taxon = ""
        if host != None and 'metadata' in hosts_dict[host]:
            print("\t".join([acc,taxon,host] + hosts_dict[host]['metadata']))
        elif host != None:
            print("\t".join([acc,taxon,host]))
    #record_data = scrape_records(gb_records)
    #taxon2family_data = fetch_taxon_family_data(record_data)
    #print(taxon2family_data)
    # for record in record_data:
    #     #print(record)
    #     if record[1] in taxon2family_data:
    #         print('\t'.join(record + taxon2family_data[record[1]]))
    #     else:
    #         print('\t'.join(record + ['', '']))

    exit()
