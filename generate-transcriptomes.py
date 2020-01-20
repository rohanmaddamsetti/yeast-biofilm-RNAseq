#!/usr/bin/env python

''' 
generate-transcriptomes.py by Rohan Maddamsetti.

use minimap2's python API (the mappy python module):
www.github.com/lh3/minimap2/python

Map the S288c reference transcriptome (from Ensembl, referenced by the kallisto docs)
to the resequenced FASTA genomes for HMY12, HMY127, HMY362.

Then use those alignments to infer transcriptomes for those three strains, using the headers
from the S288c reference transcriptome.

NOTE: some updated transcripts have bad start codons. Are these true mutations or artifacts?
      If real, then those transcripts should have little to no expression in the sleuth results.
      Do this hypothesis test as a quality control check on the genome resequencing pipeline.

'''

from os.path import dirname, realpath, join
import mappy as mp
import re

def location_from_comment(comment_str):
        ''' 
        parse the comment string to get the chromosome and location of the query.
        (rest of the header in the S288c transcriptome file).
        '''
        loc_string = comment_str.split()[1]
        _, _, chromosome, chr_start, chr_end, strand = loc_string.split(':')
        return (chromosome, chr_start, chr_end, strand)

def chromosome_from_hit(ctg):
        ''' ctg: name of the reference sequence the query is mapped to.'''
        ref_key_val, chr_key_val = ctg.split('|')
        ref = ref_key_val.split('=')[1]
        chr = chr_key_val.split('=')[1]
        return (ref,chr)


def write_dict_to_logfile(outdir, outname, my_dict, logf_suffix):
        ''' write contents of a python dictionary to a log file.'''
        logf = join(outdir,outname+logf_suffix)
        with open(logf, 'w') as fh:
                for k,v in my_dict.items():
                        fh.write(','.join([k,str(v)+'\n']))

def write_introns_to_logfile(outdir, outname, intron_dict):
        ''' write transcripts containing introns called by minimap2 to a log file.'''
        return write_dict_to_logfile(outdir, outname, intron_dict, '_intron_count.log')

def write_zero_hits_to_logfile(outdir, outname, zero_hits):
        ''' write transcripts that don't match to a log file. '''
        return write_dict_to_logfile(outdir, outname, zero_hits, '_zero_hit.log')

def write_multiple_hits_to_logPAF(a, S288c_transcripts, outdir, outname, multiple_hits):
        ''' 
        write transcripts that match to multiple times on the same chromosome to a log file
        in PAF format. Use paftools view outdir/outname.paf to visualize those alignments.
        '''
        with open(join(outdir,outname+'_multiple_hits.log.paf'),'w') as multiple_fh:
                for name, seq, _, comment in mp.fastx_read(S288c_transcripts, read_comment=True):
                        if name not in multiple_hits: continue
                        chromosome, chr_start, chr_end, strand = location_from_comment(comment)
                        for h in a.map(seq, cs=True): # traverse hits
                                _, hit_chromosome = chromosome_from_hit(h.ctg)
                                if hit_chromosome != chromosome:
                                        continue
                                outstring = name + '\t' + str(len(seq)) + '\t' + str(h) + '\n'
                                multiple_fh.write(outstring)
                                
def count_introns(h):
        ''' count the number of 'N' characters in the CIGAR string of the hit.'''
        return h.cigar_str.count('N')

def update_transcript(a, h):
        ''' 
        for each matched regex group in the cs string, concatenate to a string
        representing the aligned sequence from the reference.
        for matches, insertions, and mismatches, add to the string. 
        for deletions (in the reference), do nothing.
        skip over introns called by minimap2.
        '''

        my_transcript = ""
        idx = 0 

        ## retrieve the matching subsequence from the index.
        s = a.seq(h.ctg, h.r_st, h.r_en)

        ''' 
        see the documentation for the cs string regex at the minimap2 manpage, 
        online at: https://lh3.github.io/minimap2/minimap2.html 
        '''
        cs_regex = re.compile(r'(=[ACGTN]+|:[0-9]+|\*[acgtn][acgtn]|\-[acgtn]+|\+[acgtn]+|~[acgtn]{2}[0-9]+[acgtn]{2})')

        for m in cs_regex.findall(h.cs):
                if m.startswith(':'): ## get the length of the match and add to the transcript.
                        seq_match_len = int(m[1:])
                        my_transcript = my_transcript + s[idx:idx + seq_match_len]
                        idx = idx + seq_match_len
                elif m.startswith('*'): ## mismatch of a single base. add to the transcript.
                        my_transcript = my_transcript + s[idx:idx + 1]
                        idx = idx + 1
                elif m.startswith('-'):
                        ## deletion in the transcript; add the genome sequence to the transcript.
                        deletion_len = len(m[1:])
                        my_transcript = my_transcript + s[idx:idx + deletion_len]
                        idx = idx + deletion_len
                elif m.startswith('+'):
                        ## insertion in the transcript that is not in the genome. skip over.
                        insertion_len = len(m[1:])
                        continue
                elif m.startswith('~'): ## in an intron: skip over genome sequence and update idx.
                        ## because pattern is something like '~ag123ag' where 123 is intron length, with 4bp of flanking splice signal 'ag'.
                        intron_len = int(m[3:-2]) + 4
                        idx = idx + intron_len
                else: ## note that we don't deal with '=' for the long-form cs string.
                      raise Exception("ERROR: failed to match cs string group '{}'".format(m))

        ## check the strand, take the reverse complement if needed.
        if h.strand == -1: my_transcript = mp.revcomp(my_transcript)

        return my_transcript

        
def update_comment(comment, h, strain_name):
        ''' update the relevant fields in the S288c transcriptome annotation 
            to match the strain genome.'''
        _, hit_chromosome = chromosome_from_hit(h.ctg)
        info_string = ':'.join(['chromosome', strain_name, hit_chromosome,  str(h.r_st), str(h.r_en), str(h.strand)])
        comment_fields = comment.split()
        comment_fields[1] = info_string
        updated_comment = ' '.join(comment_fields)
        return(updated_comment)


def create_transcriptome(fasta_genome_f, S288c_transcripts, outdir, outname):
        ## do a spliced alignment.
        a = mp.Aligner(fasta_genome_f, preset='splice')
        if not a: raise Exception("ERROR: failed to load/build index file '{}'".format(fasta_genome_f))

        hit_counts = {}
        intron_counts = {}

        ## write updated transcriptome to a FASTA file. This is what we want.
        updated_transcriptome_fh = open(join(outdir, outname + '_transcriptome.fa'), 'w')
        
        ## write alignments to a PAF file. This is for internal quality control.
        ## can view with "paftools.js view outname_hits.log.paf".
        paf_fh = open(join(outdir, outname + '_hits.log.paf'), 'w')

        ## log updated transcripts with S288c query and the genome alignment for internal quality control.
        fasta_log_fh = open(join(outdir, outname + '.log.fasta'), 'w')
        
        ## log updated transcripts that don't start with an 'ATG' codon for internal quality control.
        ATG_fh = open(join(outdir, outname + '_bad_start_codons.log.fasta'), 'w')
        
        for name, seq, _, comment in mp.fastx_read(S288c_transcripts, read_comment=True): 
                chromosome, chr_start, chr_end, strand = location_from_comment(comment)
                hit_counts[name] = 0 ## initialize
                for h in a.map(seq, cs=True): # traverse hits
                        ## check that the hit is on the correct chromosome.
                        _, hit_chromosome = chromosome_from_hit(h.ctg)
                        if hit_chromosome != chromosome: continue

                        ## keep track of genes with secondary alignments, but
                        hit_counts[name] += 1
                        if h.is_primary == False: continue ## only work with primary alignments.
                        
                        ## get the updated transcript from the CIGAR or CS string of the hit
                        updated_transcript = update_transcript(a, h)
                        updated_comment = update_comment(comment, h, outname)
                        
                        official_updated_transcript_fasta = '>'+ name + ' ' + updated_comment + '\n' + updated_transcript + '\n\n'
                        updated_transcriptome_fh.write(official_updated_transcript_fasta)
                        
                        S288c_query_fasta = '>S288c '+ name + '\n' + seq + '\n'
                        updated_transcript_fasta = '>' + outname + ' ' + name + '\n' + updated_transcript + '\n'
                        genome_alignment_fasta = '>' + h.ctg + ':' + str(h.r_st) + ':' + str(h.r_en) + '\n' + a.seq(h.ctg, h.r_st, h.r_en) + '\n'

                        fasta_log_fh.write(S288c_query_fasta)
                        fasta_log_fh.write(updated_transcript_fasta)
                        fasta_log_fh.write(genome_alignment_fasta + '\n')

                        if not updated_transcript.startswith('ATG'): ## log transcripts with bad start codons.
                                ATG_fh.write(S288c_query_fasta)
                                ATG_fh.write(updated_transcript_fasta)
                                ATG_fh.write(genome_alignment_fasta + '\n')

                        PAFstring = name + '\t' + str(len(seq)) + '\t' + str(h) + '\n'
                        paf_fh.write(PAFstring)

                        ## count introns that minimap2 calls.
                        introns = count_introns(h)
                        if introns > 0: intron_counts[name] = introns

        updated_transcriptome_fh.close()
        paf_fh.close()
        fasta_log_fh.close()
        ATG_fh.close()
                        
        multiple_hits = {}
        zero_hits = {}
        for k,v in hit_counts.items():
                if v > 1:
                        multiple_hits[k] = v
                elif v == 0:
                        zero_hits[k] = v

        ## write names of transcripts that don't match to a log file.
        write_zero_hits_to_logfile(outdir, outname, zero_hits)

        ## write transcripts with introns called by minimap2 to a log file.
        write_introns_to_logfile(outdir, outname, intron_counts)
        
        ## write transcripts that match to multiple locations to a log file.
        if len(multiple_hits): write_multiple_hits_to_logPAF(a, S288c_transcripts, outdir, outname, multiple_hits)


        
def main():

    srcdir = dirname(realpath(__file__))
    assert srcdir.endswith('src')
    projdir = dirname(srcdir)
    assert projdir.endswith('yeast-RNAseq')

    datadir = join(projdir,"data")
    S288c_transcriptome_f = join(datadir,"Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa")
    HMY12_fasta_f = join(datadir,"HMY12.fasta")
    HMY127_fasta_f = join(datadir,"HMY127.fasta")
    HMY362_fasta_f = join(datadir,"HMY362.fasta")

    outdir = join(projdir,"results","transcriptomes")
    
    create_transcriptome(HMY12_fasta_f, S288c_transcriptome_f, outdir,"HMY12")
    create_transcriptome(HMY127_fasta_f, S288c_transcriptome_f, outdir,"HMY127")
    create_transcriptome(HMY362_fasta_f, S288c_transcriptome_f, outdir,"HMY362")
                        
if __name__ == "__main__":
	main()
