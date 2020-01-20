#!/usr/bin/env python

''' 
bed_to_transcriptome.py by Rohan Maddamsetti. 

This script takes yeast genome annotation in the format of a *.bed file and
a reference genome in the format of a *.fasta file
and generates a *.fasta transcriptome (set of all possible transcripts),
combinatorially generated through combinations of exons,
that will used by as input for transcript quantification with kallisto.

IMPORTANT: the next step will be to verify the quality of the *.bed annotation
and the *.fasta reference genome that are inputs to this script.

IMPORTANT: ADRIANA'S CODE HAS A FENCEPOST ERROR, WHERE THE START POSITION IS OFF BY ONE.
           FIX THIS BUG, THEN REMOVE THE COMPENSATORY HACK THAT FIXES THIS BUG.

'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from pprint import pprint

## for running and parsing pairwise alignments with MAFFT.
from Bio.Align.Applications import MafftCommandline
from io import StringIO
from Bio import AlignIO


########## HACK TO FIX FENCEPOST ERROR IN BED FILE.
def bed_to_annotation_list(bedfile):
    ''' 
    read names and positions from bed file.
    Note that positions in bed-format are zero-indexed and do not include the right
    end position, so it plays well with python out of the box.

    Input: a bed file.
    Output: a dictionary, indexed by the yeast chromosome name.
            The value of each entry is a list of tuples, where each tuple is:
            gene name, gene start, gene end (in python-friendly chromosome
            coordinates), strand.
    '''

    positions = defaultdict(list)
    with open(bedfile) as f:
        for line in f:
            line = line.strip()
            chrom, chromStart, chromEnd, name, score, strand, feattype, is_characterized  = line.split()
            chrname = chrom.rstrip('|') ## get rid of trailing '|'
            ########## SUBTRACT ONE AS A HACK TO FIX FENCEPOST ERROR IN BED FILE.
            positions[chrname].append((name, int(chromStart) - 1, int(chromEnd), strand))
    return(positions)

def fasta_to_dict(fastafile):
    ''' parse fasta file and turn into a dictionary.'''
    records = SeqIO.to_dict(SeqIO.parse(open(fastafile,'r'), 'fasta'))
    return(records)

def get_genes(bed_ref_tuple):

    ## pass in a tuple-- the user MUST make sure that the bed file
    ## and reference file correspond to each other.
    my_bed_file, my_reference = bed_ref_tuple
    
    my_bed = bed_to_annotation_list(my_bed_file)
    my_chromosomes = fasta_to_dict(my_reference)
    ''' These are key-value pairs, indexed from '1' to '17',
    where '1' is the first chromosome, '16' is the last
    chromosome, and '17' is the mitochondrion. '''

    ''' my_bed is keyed to the name of the Genbank reference for each
    S. cerevisiae S288c chromosome. By contrast, my_chromosome
    is indexed by the chromosome number from '1' to '17' where
    '17' is the mitochondrion.
    '''
    chrindex_dict = {x:str(i+1) for i,x in enumerate(my_bed.keys())}
        
    ## extract genes from the chromosome.
    short_seq_record_dict = {}
    for chrname in my_bed:
        for (gene_name, start, stop, strand) in my_bed[chrname]:
            long_seq_record = my_chromosomes[chrindex_dict[chrname]]
            long_seq = long_seq_record.seq
            alphabet = long_seq.alphabet
            short_seq = Seq(str(long_seq)[start:stop], alphabet)
            if strand == '-': ## then take the reverse complement.
                short_seq = short_seq.reverse_complement()
            short_seq_record = SeqRecord(short_seq, id=gene_name, description='')
            short_seq_record_dict[gene_name] = short_seq_record
    return(short_seq_record_dict)

def align_strain_gene_to_S288c_transcript(gene_dict,ref_transcript_dict, bedf):
    for k,v in gene_dict.items():
        if k not in ref_transcript_dict:
            print('WARNING: gene', k, 'in bed file', bedf, 'not found in S288')
        else: ## align strain sequence against the S288c transcript.
            strain_gene = gene_dict[k]
            ref_transcript = ref_transcript_dict[k]

            ## write genes to align to file.
            ## IMPORTANT: THERE ARE SERIOUS DISCREPANCIES BETWEEN MANY OF THESE SEQUENCES--
            ## THIS LOOKS LIKE A BUG IN ADRIANA'S ANNOTATION CODE TO ME...
            aln_me_file = "../temp/" + k + ".fasta"
            with open(aln_me_file, "w") as output_handle:
                SeqIO.write([ref_transcript, strain_gene], output_handle, "fasta")
            mafft_cline = MafftCommandline(input=aln_me_file)
            ## run mafft and parse the input.
            stdout, stderr = mafft_cline()
            align = AlignIO.read(StringIO(stdout), "fasta")

            ## print alignment for debugging purposes.
            ## I am debugging as follows: python bed_to_transcriptome.py > ../test_alignments.txt
            print("Alignment length %i" % align.get_alignment_length())
            for record in align:
                print("%s - %s" % (record.seq, record.id))
                print()
            print('**************************')
            ##quit()

def main():

    ## load the S288c reference transcriptome,
    ## and store as dict of gene_name : SeqRecord pairs.
    ## the S288c transcriptome is from Ensembl, via the kallisto website.
    S288c_transcriptome_f = '../data/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa'
    ref_transcript_dict = {k.split('_mRNA')[0]:v for k,v in fasta_to_dict(S288c_transcriptome_f).items()}
    
    ## load the resequenced genomes for HMY12, HMY127, HMY362.
    input_tuples = [('../data/HMY12NewBedFile_OrigChrom.bed',
                     '../data/HMY12_NewReference.fasta'),
                    ('../data/HMY127NewBedFile_OrigChrom.bed',
                     '../data/HMY127Genome.fasta'),
                    ('../data/HMY362NewBedFile_OrigChrom.bed',
                     '../data/HMY362.genome.fasta')]

    for tup in input_tuples:
        assert tup[0].endswith('.bed') and tup[1].endswith('.fasta')
        ## transcribe strain genes and store as dict of gene_name : SeqRecord pairs.
        gene_dict = get_genes(tup)

        ## pass bed file into this function to report transcripts
        ## that are missing from S288c reference transcriptome.
        align_strain_gene_to_S288c_transcript(gene_dict, ref_transcript_dict, tup[0])                
        
        quit() ## There are inconsistencies across these inputs from Adriana!!!
        ## before running on everything, re-run the code that generated the *.bed and *.fasta
        ## files in the first place to fix those inconsistencies.

        ## right now, this code only works for HMY12.

main()

###################
## extra code to put back in later follows.
quit()

# write to file
with open('output.fasta', 'w') as f:
    SeqIO.write(short_seq_records, f, 'fasta')
