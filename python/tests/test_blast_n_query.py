'''
Created on Jan. 30, 2023

@author: pleyte
'''
import unittest
from Bio.Blast import NCBIXML
from reproduce.pmid33033132 import blast_n_query


class Test(unittest.TestCase):


    def test_extract_taxa(self):
        query_results = NCBIXML.parse(open('../../galaxy/test-data/SRR12352293_1.blastn.xml'))
        taxa = blast_n_query._extract_taxa(query_results, float('1e-20'), 'xyz')    
        query_results.close()
        self.assertEqual(193, len(taxa), "Taxa list size")

    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()