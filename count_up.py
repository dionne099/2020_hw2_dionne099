#!/usr/bin/env python3

# this is a python script template
# used curl -O in bash to obtain data files

gff="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz"
fasta="Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz"

import os,gzip,itertools,csv,re,sys

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

def my_custom_row_func(row): 
    start = ''
    end = ''

def aspairs(f):
    seq_id = ''
    sequence = ''
    for header,group in itertools.groupby(f, isheader):
        if header:
            line = next(group)
            seq_id = line[1:].split()[0]
        else:
            sequence = ''.join(line.strip() for line in group)
            yield seq_id, sequence
            
if not os.path.exists(gff):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/gff3/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.37.gff3.gz")

if not os.path.exists(fasta):
    os.system("curl -O ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/dna/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.dna.chromosome.Chromosome.fa.gz")

print("/////////////////////////")
            
count_up_genes = 0
total_length = 0
features = 0 

with gzip.open(gff,"rt") as fh: 
    gff = csv.reader(fh,delimiter="\t") 
    for row in gff: 
        if row[0].startswith("#"): 
            continue 
        if "gene" == row[2]: 
            count_up_genes += 1                                    
            Exonlen = int(row[4]) - int(row[3]) #this works!
            total_length += Exonlen 
            
            
print("There are",count_up_genes ,"genes in the GFF file.")                     
print("The total length of the genes in the GFF file are",total_length, "bp.")

print("===========================")

f_len=0

with gzip.open(fasta, "rt") as fh:
    
    seqs = aspairs(fh)
    
    for seq in seqs:
        f_len = len(seq[1])
        
print("The length of the chromosome is",f_len,"bp.")   
print("===========================")

cds= 100*(total_length/f_len)
print("The % coding is", cds, "%")
