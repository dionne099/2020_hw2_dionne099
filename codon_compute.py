#!/usr/bin/env python3

import os, gzip, itertools 
from collections import Counter 

# this is code which will parse FASTA files
# define what a header looks like in FASTA format
def isheader(line):
    return line[0] == '>'

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
    
url1="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2/cds/Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
url2="ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_0_collection/mycobacterium_tuberculosis_h37rv/cds/Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"
file1="Salmonella_enterica_subsp_enterica_serovar_typhimurium_str_lt2.ASM694v2.cds.all.fa.gz"
file2="Mycobacterium_tuberculosis_h37rv.ASM19595v2.cds.all.fa.gz"

if not os.path.exists(file1):
    os.system("curl -O -L %s"%(url1))

if not os.path.exists(file2):
    os.system("curl -O -L %s"%(url2))
    
#find the total number of genes and total length for both files 

print("///////////////////")

gene_num_1= 0 
gene_length_1= 0
total_gene_length= 0
dna=[]

with gzip.open(file1,"rt") as fh:
    
    seqs = aspairs(fh)
    #print(seqs)
    
    for seq in seqs:
        gene_num_1 += 1
        dna +=list(seq[1])
        gene_length_1 = len(seq[1])
        total_gene_length += gene_length_1
        #print(seqs.keys())
        
print("There are", gene_num_1, "genes in file 1")
print("The total length of these genes are", total_gene_length, "bp.")
#print(len(dna))
print("===================")

gene_num_2= 0 
gene_length_2= 0
total_gene_length_2= 0
dna2= []

with gzip.open(file2,"rt") as fh:
    seqs = aspairs(fh)
    
    for seq in seqs:
        gene_num_2 += 1 
        dna2 +=list(seq[1])
        gene_length_2 = len(seq[1])
        total_gene_length_2 += gene_length_2
        
print("There are", gene_num_2, "genes in file 2")
print("The total length of these genes are", total_gene_length_2, "bp.")
#print(len(dna2))
print("===================")

#find the total length of both files 

combined_gene_length= total_gene_length + total_gene_length_2
print("The combined gene lengths of both files are:",combined_gene_length,"bp.")

#create and combine the values of the dictionaries to get total G/ C content

bases = {'A':0, 'C':0, 'G':0, 'T':0 }
for l in dna:
   bases[l] += 1

#print(bases)

bases2 = {'A':0, 'C':0, 'G':0, 'T':0 }
for l in dna2:
   bases2[l] += 1

#print(bases2)

total_bases =  {x: bases.get(x, 0) + bases2.get(x, 0) 
                    for x in set(bases).union(bases2)} 
#print("The dictionary for both files:", total_bases)
#print("===================")

GCKeys = ["G", "C"]
ATKeys = ["A", "T"]

filterByKey = lambda keys: {x: total_bases[x] for x in keys}
GCdict= filterByKey(GCKeys)
ATdict= filterByKey(ATKeys)

GCcontent = GCdict.values()
ATcontent = ATdict.values()

GCtotal = sum(GCcontent)
ATtotal = sum(ATcontent)
GCATtotal = GCtotal + ATtotal

GCpercent = ((GCtotal/GCATtotal)*100)
print("The GC Content is:", GCpercent, "%")
print("===================")

totalnumcodons_gen1= len(dna)/3
totalnumcodons_gen2= len(dna2)/3
print("The total number of codons in genome 1 is", totalnumcodons_gen1)
print("The total number of codons in genome 2 is", totalnumcodons_gen2)
print("===================")

#mission create a codon table 

codon1 = {}
codon_list1= [] 
i=0

with gzip.open(file1,"rt") as fh:
    seqs = dict(aspairs(fh))
    for seqname in seqs:
        #print(seqs[seqname])
        #if i == 0:
            #print("seq is {}".format(seq))
            #i += 1
            for i in range(0,len(seqs[seqname]),3):
                c = seqs[seqname][i:i+3]
                if c in codon1:
                    codon1[c] += 1 #at whatever position c is, in codon, add 1 to the value  
                else:
                    codon1[c] = 1 #at whatever position c is, set it to 1
            codon_list1 = [seqs[seqname][i:i+3] for i in range(0, len(seqs[seqname]), 3)] 
for c in codon_list1:
             # now 
    if c in codon1:
        codon1[c] += 1
    else:
        codon1[c] = 1
                        
#print(codon1)

codon2 = {}
codon_list2= [] 
i=0

with gzip.open(file2,"rt") as fh:
    seqs = dict(aspairs(fh))
    for seqname in seqs:
        #print(seqs[seqname])
        #if i == 0:
            #print("seq is {}".format(seq))
            #i += 1
            for i in range(0,len(seqs[seqname]),3):
                c = seqs[seqname][i:i+3]
                if c in codon2:
                    codon2[c] += 1 #at whatever position c is, in codon, add 1 to the value  
                else:
                    codon2[c] = 1 #at whatever position c is, set it to 1
            codon_list2 = [seqs[seqname][i:i+3] for i in range(0, len(seqs[seqname]), 3)] 
for c in codon_list2:
             # now 
    if c in codon2:
        codon2[c] += 1
    else:
        codon2[c] = 1
                        
#print(codon2)

import pandas as pd 

df1 = pd.DataFrame(list(codon1.items()), columns= ['Codon','Species 1 Frequency'])
df2 = pd.DataFrame(list(codon2.items()), columns= ['Codon','Species 2 Frequency'])
#print(df1)
#print(df2)

df3 = pd.merge(df1, df2)
print(df3)
