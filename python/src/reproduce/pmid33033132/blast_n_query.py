#!/usr/local/bin/python3
# encoding: utf-8
'''
This script takes a fasta file, queries blastn, and writes out the matching species associated with each sequence. 

@author:     jpleyte
@copyright:  2022 Jay Pleyte. All rights reserved.

@license:    GNU GENERAL PUBLIC LICENSE v3

@deffield    updated: 2022.11.13
'''

import argparse
import logging
import sys
import time
import csv
import re

from Bio.Blast import NCBIWWW, NCBIXML
from reproduce.pmid33033132.taxon import Taxon
from Bio import Entrez, SeqIO

__version__ = "0.0.4"

logger = logging.getLogger()
streamHandler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
streamHandler.setFormatter(formatter)
logger.addHandler(streamHandler)
logger.setLevel(logging.DEBUG)

def _query(fasta_file:str):
    with open(fasta_file, mode="r") as f:
        sequence_data = f.read()
        
    logger.debug("Querying blastn server")
    start_time = time.perf_counter()
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)        
    end_time = time.perf_counter()
    num_reads = sequence_data.count('\n')
    logger.debug(f"blastn query containing {num_reads} sequences completed in {end_time - start_time:0.2f} seconds")

    results = NCBIXML.parse(result_handle)

    # For testing, uncomment the following code to write the results to an xml file 
    # results = None
    # blast_results = result_handle.read()
    # logger.debug(f"Saving results to xml file")     
    # with open('/tmp/result.xml', 'w') as save_file:
    #     save_file.write(blast_results)
    
    return results

def _parse_species(blastn_hitdef: str):
    '''
    Parse the species out of a blastn Hit_def string 
    '''
    
    # Barbus barbus genome assembly, chromosome: 24
    match = re.search("^(.*) genome assembly, chromosome:.*$", blastn_hitdef)    
    if match:
        return match.group(1).strip()
    
    # Example: PREDICTED: Polypterus senegalus lysosomal thioesterase PPT2-A-like (LOC120539052), transcript variant X5, mRNA
    match = re.search("^PREDICTED:(\b.*)uncharacterized.*", blastn_hitdef)
    if match:
        return match.group(1).strip()
    
    # Example: Arcobacter venerupis strain LMG 26156 chromosome, complete genome
    match = re.search("^(.*) strain .* chromosome.*", blastn_hitdef)
    if match:
        return match.group(1).strip()
    
    # Example: Shewanella baltica strain NCTC10737 genome assembly, chromosome: 1
    match = re.search("^(.*) strain .* genome assembly, chromosome.*", blastn_hitdef)
    if match:
        return match.group(1).strip()
    
    # Example: Massilia sp. NP310 chromosome, complete genome
    match = re.search("^(.*) chromosome, complete genome", blastn_hitdef)
    if match:
        return match.group(1).strip()
    
    # Example: Shewanella baltica OS678, complete genome
    match = re.search("^(.*), complete genome", blastn_hitdef)
    if match:
        return match.group(1).strip()
    
    # Example: Gossypioides kirkii chromosome KI_11
    match = re.search("^(.*) chromosome .*$", blastn_hitdef)
    if match:
        return match.group(1).strip()
    
    # Example: Oppiella nova
    match = re.search("^([A-Z][a-z]+\s[a-z]+)$", blastn_hitdef)
    if match:
        return match.group(1).strip()
    
    logger.warning(f"Unable to parse species from string: {blastn_hitdef}")
    return None
    
def _extract_taxa(blastn_records, threshold:float, sample: str):
    '''
    Iterate over the blastn results, filter by threshold, and return a list of Taxon object. 
    ''' 
    taxa = []
    
    n = 0
    n_kept = 0
    
    for record in blastn_records:
        if record.alignments: 
            for align in record.alignments: 
                for hsp in align.hsps:
                    n += 1
                    if hsp.expect < threshold:
                        n_kept += 1
                        taxon = Taxon()                        
                        taxon.accession = align.accession
                        taxon.e_value = hsp.expect
                        taxon.hitdef = align.hit_def
                        taxon.sample = sample
                        taxon.species = _parse_species(taxon.hitdef)
                        taxa.append(taxon)
    
    logger.info(f"Keeping {n_kept} of {n} matches that meet threshold of {threshold}.")

    return taxa

def _lookup_domain(tax_csv, species):
    '''
    Search the taxid2parents.tsv file for the domain (aka superkingdom) of a species. 
    '''
    for row in tax_csv:
        if len(row) < 4:
            continue
        
        # row[0] is the id, and row[1] is the species
        (row_species, label, species_id) = row[1].split('|')
        
        # Sometimes the row has a strain identifier before the species
        row_strain = None
        if label  == 'strain':
            row_strain = row_species
            (row_species, label, species_id) = row[2].split('|')
        
        if(label != 'species'):
            continue

        if species == row_species or row_strain == species:
            # domain is always the second to last column 
            (domain, label, domain_id) = row[len(row)-2].split('|')
            assert(label == 'superkingdom')
            return domain
    
    logger.warning(f"Unable to find species in mapping file: {species}")

def _update_domain(taxid2parents_file: str, taxa: list):
    '''
    Find species in the taxid2parents file in order to determine the domain
    ''' 
    species_to_domain_map = {}
    
    with open(taxid2parents_file, mode="r") as tax_file:
        tax_csv = csv.reader(tax_file, delimiter="\t")

        line_n = 0
        for taxon in taxa:
            line_n = line_n + 1
            if taxon.species in species_to_domain_map:
                taxon.domain = species_to_domain_map[taxon.species]
            elif taxon.species:
                taxon.domain = _lookup_domain(tax_csv, taxon.species)
                species_to_domain_map[taxon.species] = taxon.domain
                tax_file.seek(0)
            
            if line_n % 10 == 0:
                logger.debug(f"Updating domain in taxon {line_n} of {len(taxa)}")
        

        
def _write(out_file: str, taxa: list):
    '''
    '''
    fields = ['sample', 'accession', 'e_value', 'hitdef', 'species', 'domain']
    with open(out_file, 'w') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(fields)
        logger.debug(f"Writing taxa to {out_file}")
        for taxon in taxa:
            csv_writer.writerow([taxon.sample, taxon.accession, taxon.e_value, taxon.hitdef, taxon.species, taxon.domain])

def _parse_args():    
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='A test of how scripts expect individual file parameters handle Galaxy datasets and collections.')
    
    parser.add_argument('-i', '--in_file',  
                        help='Input FASTA file',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('-o', '--out_file',  
                        help='Output csv file',
                        type=argparse.FileType('w'),
                        required=True)
    
    parser.add_argument('-t', '--threshold',  
                        help='blastn e-value maximum',
                        required=False,
                        default="1e-20")

    parser.add_argument('-s', '--sample',
                        help='sample name',
                        required=False,
                        default="")
    
    parser.add_argument('-x', '--taxid2parents',
                        help='TSV file containing mappings from species to domain (see https://www.genome.jp/tools-bin/taxsummary)"',
                        type=argparse.FileType('r'),
                        required=True)
    
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    
    args = parser.parse_args()
    return args
                        
def _main():
    '''
    main function
    '''
    args = _parse_args()
            
    query_results = _query(args.in_file.name)
    
    # Extract accession id from query results
    taxa = _extract_taxa(query_results, float(args.threshold), args.sample)
    query_results.close()
    
    # Update the taxa objects with domain
    _update_domain(args.taxid2parents.name, taxa)
    
    _write(args.out_file.name, taxa)
    
if __name__ == "__main__":
    _main()