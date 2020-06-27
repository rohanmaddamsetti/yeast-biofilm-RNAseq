#!/usr/bin/env python

'''
compare-GSNAP-kallisto-counts.py by Rohan Maddamsetti.

check the raw read counts in results/GSNAP-quantification-comparison/Biofilm.counts.txt
against the transcript counts made by kallisto.

prints out tables for analysis with compare-GSNAP-kallisto-counts.R.

also-- check whether rRNA accounts for the mapping discrepancy.
To do so, compare the Gene IDS in Biofilm.counts.txt against the
transcriptome that I used for kallisto.
'''

import os
import itertools
import subprocess

def compare_mapped_genes():
    ## slurp in a list of Gene IDs in Biofilm counts.
    fh1 = open("../results/GSNAP-quantification-comparison/Biofilm.counts.txt", "r")

    GSNAP_ids = []
    for i,l in enumerate(fh1):
        if i == 0: continue
        l.strip()
        ldata = l.split()
        if len(l) > 1:
            GSNAP_ids.append(ldata[0])
    fh1.close()

    fh2 = open("../data/Saccharomyces_cerevisiae.R64-1-1.cdna.all_AND_FLO8.fa")

    transcript_ids = []
    for j,l in enumerate(fh2):
        if not l.startswith('>'): continue
        ldata = l.split()
        gene_field = ldata[3]
        gene = gene_field.split(':')[-1]
        transcript_ids.append(gene)

    s1 = set(GSNAP_ids)
    s2 = set(transcript_ids)

    print("number of genes in GSNAP but not transcriptome:", len(s1-s2))
    print("number of genes in transcriptome but not GSNAP:", len(s2-s1))


def cross_check_read_counts():
    ## let's count up the mapped reads in Biofilm.counts.txt, and cross-check
    ## against the mapped numbers.
    fh = open("../results/GSNAP-quantification-comparison/Biofilm.counts.txt", "r")

    read_counts = {}
    idx_to_sample = {}
    for i,l in enumerate(fh):
        l.strip()
        ldata = l.split()
        if i == 0: ## map index of each field to sample name.
            idx_to_sample = {j:x for j,x in enumerate(ldata)}
            read_counts = {x:0 for x in ldata[1:]} ## skip gene_id field
            continue
        if len(l) > 1: ## skip empty lines due to weird newline char
            for j,datum in enumerate(ldata):
                if j == 0: continue ## skip gene_id field
                sample_name = idx_to_sample[j]
                read_counts[sample_name] = read_counts[sample_name] + int(datum)
    print(read_counts)


def write_raw_count_table():
    ## let's write out a nicely formatted table to compare the quantifications
    ## made by GSNAP and kallisto.

    outfh = open("../results/GSNAP-quantification-comparison/GSNAP_raw_counts.csv", "w")
    
    ## print header for this table.
    header = "gene_id,sample,raw_counts\n"
    outfh.write(header)
    
    fh3 = open("../results/GSNAP-quantification-comparison/Biofilm.counts.txt", "r")

    idx_to_sample = {}
    for i,l in enumerate(fh3):
        l.strip()
        ldata = l.split()
        if i == 0: ## map index of each field to sample name.
            idx_to_sample = {j:x for j,x in enumerate(ldata)}
            continue
        if len(l) > 1: ## skip empty lines due to weird newline char
            for j, datum in enumerate(ldata):
                if j == 0:
                    ## save the geneid and move on to the other fields.
                    geneid = datum
                    continue 
                sample_name = idx_to_sample[j]
                ## print the row.
                outfh.write(','.join([geneid,sample_name,datum]))
                outfh.write("\n")
    outfh.close()

def cp_kallisto_abundance():
    ''' all this does is copy the relevant kallisto abundance data into 
     results/GSNAP-quantification-comparison. 
    pandas is too much of a pain in the ass to use here. '''

    strains = ['HMY12','HMY127','HMY362']
    for s in strains:
        ## get the kallisto results for this strain, using this reference genome.
        kallisto_ref_dir = s + '-ref-kallisto-output'
        samples = [s + ''.join(m) for m in itertools.product(('a','b','c'), ('LD','YPD'))]
        for sam in samples:
             my_abundance_cp_in = os.path.join('../results',kallisto_ref_dir,sam,'abundance.tsv')
             my_outf = sam + '-abundance.tsv'
             my_abundance_cp_out = os.path.join('../results','GSNAP-quantification-comparison',my_outf)
             subprocess.run(['cp',my_abundance_cp_in, my_abundance_cp_out])
             
            
        
                
def main():
    write_raw_count_table()
    cp_kallisto_abundance()

main()
    
    
    
