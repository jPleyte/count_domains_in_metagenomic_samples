#!/usr/local/bin/python3
# encoding: utf-8
'''
This script takes a fasta file and writes out the size of the file

@author:     jpleyte
@copyright:  2022 Jay Pleyte. All rights reserved.

@license:    GNU GENERAL PUBLIC LICENSE v3

@contact:    jpleyte@users.noreply.github.com
@deffield    updated: 2022.11.13
'''

import argparse
import logging
import sys
import os

__version__ = "0.0.2"

logger = logging.getLogger()
streamHandler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
streamHandler.setFormatter(formatter)
logger.addHandler(streamHandler)
logger.setLevel(logging.DEBUG)

def _process(in_file:str, out_file:str):
    file_size = os.path.getsize(in_file)
    logger.debug(f"{in_file} is of size {file_size}. Writing to {out_file}")
    
    with open(out_file, mode="w") as f:
        f.write(f"{in_file} is of size {file_size}\n")

def _parse_args():    
    '''
    Validate and return command line arguments.
    '''
    parser = argparse.ArgumentParser(description='A test of how scripts expect individual file parameters handle Galaxy datasets and collections.')
    
    parser.add_argument('-i', '--in_file',  
                        help='Input FASTA',
                        type=argparse.FileType('r'),
                        required=True)

    parser.add_argument('-o', '--out_file',  
                        help='Output text file',
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