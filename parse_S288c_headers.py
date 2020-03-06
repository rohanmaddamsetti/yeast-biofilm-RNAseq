#!/usr/bin/env python

''' 
parse_S288c_headers.py by Rohan Maddamsetti.

This script parses the headers for the S288c reference transcriptome in
data/Saccharomyces_cerevisiae.R64-1-1.cdna.all_AND_FLO8.fa .

print out a "csv" file, using the pipe character '|' as the separator
character. This is because the description field has commas and other
punctuation.

Usage: python parse_S288c_headers.py > ../results/transcript_annotation.csv

'''

import re

def main():

    f = "../data/Saccharomyces_cerevisiae.R64-1-1.cdna.all_AND_FLO8.fa"

    regex = re.compile("^\>([\w-]+) cdna chromosome:R64-1-1:([\w:-]+) gene:([\w-]+) gene_biotype:(\w+) transcript_biotype:(\w+) (.*)description:(.+)$")
    
    
    with open(f,"r") as fh:
        print("mRNA|chromosome|chr_start|chr_end|strand|gene|gene_biotype|transcript_biotype|gene_symbol|description")
        for l in fh:
            if not (l.startswith('>')): continue
            l = l.strip()
            m = regex.match(l)
            if m:
                mRNA = m.group(1)
                chromosome_location = m.group(2)
                chromosome, chr_start, chr_end, strand = chromosome_location.split(':')
                gene = m.group(3)
                gene_biotype = m.group(4)
                transcript_biotype = m.group(5)
                gene_sym_str = m.group(6)
                description = m.group(7)
                gene_symbol = 'NA'
                if len(gene_sym_str): gene_symbol = gene_sym_str.split(':')[-1]
                assert '|' not in description
                print('|'.join([mRNA,chromosome,chr_start,chr_end,strand,gene,gene_biotype,transcript_biotype,gene_symbol,description]))
            else:
                print("Header not matched",l)
                quit()

main()
