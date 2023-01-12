#!/usr/local/bin/python3
# encoding: utf-8
'''
This script takes a fasta file and writes out the size of the file

@author:     jpleyte
@copyright:  2022 Jay Pleyte. All rights reserved.

@license:    GNU GENERAL PUBLIC LICENSE v3

@deffield    updated: 2022.11.13
'''

import argparse
import logging
import os
import sys
import time

from Bio.Blast import NCBIWWW

__version__ = "0.0.2"

logger = logging.getLogger()
streamHandler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
streamHandler.setFormatter(formatter)
logger.addHandler(streamHandler)
logger.setLevel(logging.DEBUG)

def _process(fasta_file:str, blastn_result_file:str):
    with open(fasta_file, mode="r") as f:
        sequence_data = f.read()
        
        logger.debug("Querying blastn server")
        start_time = time.perf_counter()
        result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)
        end_time = time.perf_counter()
        num_reads = sequence_data.count('\n')
        logger.debug(f"blastn query containing {num_reads} sequences completed in {end_time - start_time:0.2f} seconds")
        
        blast_results = result_handle.read()
                
        logger.debug(f"Saving results to {blastn_result_file}")     
        with open(blastn_result_file, 'w') as save_file:
            save_file.write(blast_results)
        
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
                        help='Output xml file',
                        type=argparse.FileType('w'),
                        required=True)

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    
    args = parser.parse_args()
    return args
    
def _main():
    '''
    main function
    '''
    args = _parse_args()
    _process(args.in_file.name, args.out_file.name)

if __name__ == "__main__":
    _main()