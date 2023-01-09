# Reproduce Results from Metagenomics of Wastewater Influent From Southern California

The goal of this project is to reproduce the results documented in the 2020 paper titled Metagenomics of Wastewater Influent from Southern California Wastewater Treatment Facilities in the Era of COVID-19 (1). In this project 

This project uses the same dataset from the Potham et al paper but uses its own set of tools to perform the analysis. For example, rather than using  Kraken 2 (2) to assign taxonomy to reads, I use Python to submit Fasta data to BLAST (3) using Biopython's Bio.Blast.NCBIWWW module (4). Each step is carried out using a Galaxy (5) pipeline; from downloading the publicly available data and performing QC, to analysing the data and generating reports. 

## Pipeline Steps

1. Download 17 archives containing fasta files

	- See https://ncbiinsights.ncbi.nlm.nih.gov/2021/01/05/bioproject-datasets/ for one way to download the sample data. 

	- See [Uploading from short read archive](https://galaxyproject.org/tutorials/collections/#uploading-from-short-read-archive) for a video on how to create a Galaxy dataset collection from data deposited to the Short Read Archive (SRA) at NCBI
		- Go to https://www.ncbi.nlm.nih.gov/sra/?term=SRP274310 | Send To | File | Accession List, and download accession list ``SraAccList.txt``

	- Consider using the [SRA Toolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) to download the fasta files: ``fastq-dump --stdout -X 16 SRR12352293``	
		- You'll want to use --split-3
		- astq-dump is  deprecated. Use fasterq-dump instead
		- ``fasterq-dump --split-3


	- If you want to compress the files, gzip is not recommend. Genozip,[14] DSRC and DSRC2, FQC, LFQC, Fqzcomp, and Slimfastq. 
	
	
2. Download paired-end fastq sequences from shotgun samples
	- Use the "Faster Download and Extract Reads in FASTQ" Galaxy tool.
3. Run the fastq files through QC (eg FastaQC)
4. If there is a way to remove bad data, do that
4. Convert fastq to fasta
2. For each Fasta file, submit each sequence to BLAST to determine the sample taxonomy
	- See this [BLAST Guide](https://guides.nnlm.gov/tutorial/ncbi-blast-finding-and-comparing-sequences/single-page) on interpretting blast. 
2. Calculate Phred quality score
3. Store each taxonomic classification in a database so they can be counted later 
4. Generate a report showing the abundances of taxa in each sample. 

## Data

17 samples of "shotgun metagenomic sequencing of wastewater"
* [SRR12352293](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR12352293&display=download)
* [SRR12352294](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR12352294&display=download)
* [SRR12352295](https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&acc=SRR12352295&display=download)
* SRR12352296, 
SRR12352297, 
SRR12352298, 
SRR12352299, 
SRR12352300, 
SRR12352301, 
SRR12352302, 
SRR12352303, 
SRR12352304, 
SRR12352305, 
SRR12352306, 
SRR12352307, 
SRR12352308, 
SRR12352309

Note: NCBI's site does this sort of analysis for you. See the "Analysis" tab: https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page_size=10&acc=SRR12352293&display=analysis
					

## Galaxy Tools

* sra_tools from  https://toolshed.g2.bx.psu.edu/ 

## References

1. Rothman JA, Loveless TB, Griffith ML, Steele JA, Griffith JF, Whiteson KL. Metagenomics of Wastewater Influent from Southern California Wastewater Treatment Facilities in the Era of COVID-19. Microbiol Resour Announc. 2020 Oct 8;9(41):e00907-20. doi: 10.1128/MRA.00907-20. PMID: 33033132; PMCID: PMC7545286.

2. Wood DE, Lu J, Langmead B. 2019. Improved metagenomic analysis with Kraken 2. Genome Biol 20:257. doi: 10.1186/s13059-019-1891-0.

3. BLAST

4. https://www.tutorialspoint.com/biopython/biopython_overview_of_blast.htm

5. Galaxy 