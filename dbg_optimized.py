##########################################################################################
#                                                                                        #
#  Python 3.6 de-Bruijn graph generator and differential expression encoder              #
#                                                                                        #
#  Base code taken from https://raw.githubusercontent.com/pmelsted/dbg/master/dbg.py     #
#                                                                                        #
#  Developed at Megaprobe-lab, PI Dr. Humberto Ortiz-Zuazaga                             #
#  Developers: Kevin Legarreta, Louis Gil, Israel O. Dilan-Pantojas                      #
#  Program description: Given Fastq files encoding sequences for two organisms, this     #
#                       program generates a de-Bruijn graph which encodes information    #
#                       about the coverage of k-sized regions across both organisms.     #
#                       This information is useful for the analysis of differential      #
#                       expression, theoretically the graph structure approach should    #
#                       recover more information about the organisms genetic make-up     #
#                       when compared to traditionally generated assembled contigs.      #        
#                                                                                        #
##########################################################################################

# Future Works ;)
# import multiprocessing as mp

#########################
#                       #
#  Module Imports       #
#                       #
#########################

import argparse
import collections
from Bio import Seq, SeqIO, SeqRecord
import gfapy

#########################
#                       #
#  Sequence Methods     #
#                       #
#########################

#  Create reverse complement
def twin(km):
    return Seq.reverse_complement(km)

# Chunk reads into k sized chunks 
def kmers(seq,k):
    for i in range(len(seq)-k+1):
        yield seq[i:i+k]

# Move sequence analysis one nucleotide forward  
def fw(km):
    for x in 'ACGT':
        yield km[1:]+x

# Move sequence analysis one nucleotide backward
def bw(km):
    for x in 'ACGT':
        yield x + km[:-1]

#########################
#                       #
#  File Read Methods    #
#                       #
#########################

# Count kmers and build dictionary with counts
# Use limit as cut off for very down regulated seq 
def build(fn,k,limit):
    d = collections.defaultdict(int)
    # For Each filename in input parse fast.q file 
    for f in fn:
        reads = SeqIO.parse(f,'fastq')
        # Extract reads 
        for read in reads:
        # Maybe it's possible to reduce memory overhead here!    
            for seq in str(read.seq).split("N"):
                # Extract kmers from forward
                for km in kmers(seq,k):
                    d[km] +=1
                # Extract kmers form reverse
                for km in kmers(twin(seq),k):
                    d[km] += 1
    # Remove k-mer elements of dict that appear less than limit times
    for i in [x for x in d if d[x] <= limit]:
            del d[i]
    return d

#########################
#                       #
#  Sequence Methods     #
#                       #
#########################    

# Merge dictionaries to create homunculus 
def merge_dicts(d1,d2):
    return {i: (d1[i], d2[i]) for i in set(list(d1.keys()) + list(d2.keys()))}

# Stringify
def contig_to_string(c):
    return (c[0] + ''.join(x[-1] for x in c[1:]))


#########################
#                       #
#  Graph Methods        #
#                       #
#########################

# Get next k in contig
def get_contig(d,km):
    # print("here")
    # return ((contig_to_string(get_contig_forward(d,km)), get_contig_forward(d,km)) if km in fw(get_contig_forward(d,km)[-1]) else ( (contig_to_string([twin(x) for x in get_contig_forward(d,twin(km))[-1:0:-1]] + get_contig_forward(d,km))), ([twin(x) for x in get_contig_forward(d,twin(km))[-1:0:-1]] + get_contig_forward(d,km))) )
    c_fw = get_contig_forward(d,km)
    c_bw = get_contig_forward(d,twin(km))
    if km in fw(c_fw[-1]):
        c = c_fw
    else:
        c = [twin(x) for x in c_bw[-1:0:-1]] + c_fw
    return contig_to_string(c),c
        
# Move from k to k+1 and break if loop encountered
def get_contig_forward(d,km):
    c_fw = [km]
    
    while True:
        if sum(x in d for x in fw(c_fw[-1])) != 1:
            break
        
        cand = [x for x in fw(c_fw[-1]) if x in d][0]
        if cand == km or cand == twin(km):
            break # break out of cycles or mobius contigs
        if cand == twin(c_fw[-1]):
            break # break out of hairpins
        
        if sum(x in d for x in bw(cand)) != 1:
            break

        c_fw.append(cand)
    # Return contig iterated forward     
    return c_fw

# Create Graph from all contigs 
def all_contigs(d,k):
    done = set()
    r = []
    for x in d:
        if x not in done:
            s,c = get_contig(d,x)
            for y in c:
                done.add(y)
                done.add(twin(y))
            r.append(s)
    
    G = {}
    heads = {}
    tails = {}
    for i,x in enumerate(r):
        G[i] = ([],[])
        heads[x[:k]] = (i,'+')
        tails[twin(x[-k:])] = (i,'-')
    
    for i in G:
        x = r[i]
        for y in fw(x[-k:]):
            if y in heads:
                G[i][0].append(heads[y])
            if y in tails:
                G[i][0].append(tails[y])
        for z in fw(twin(x[:k])):
            if z in heads:
                G[i][1].append(heads[z])
            if z in tails:
                G[i][1].append(tails[z])

    # Returns Graph and list of reads 
    return G,r

#########################
#                       #
#  Output Methods       #
#                       #
#########################

# Write graph links to file
def get_links(cs,d,k,s,fd):
    for x in range(0,len(cs)-k):
        keyA = cs[x:x+k]
        keyB = cs[(x+1):(x+1)+k]
        
        kmerA = (('%s:%s:(A:%s,B:%s)')%(s,x,d[keyA][0],d[keyA][1]))
        kmerB = (('%s:%s:(A:%s,B:%s)')%(s,x+1,d[keyB][0],d[keyB][1]))

        fd.write("L\t%s\t+\t%s\t+\t%sM\n"%(kmerA,kmerB,(k-1)))

# Write kmers to file
def get_kmers(cs,d, k, s,fd):
    global g, listofkmers,listoflinks,lastkmerid

    for x in range(0,len(cs)+1-k):              #  to get all subsegment, holds the all the kmers
        key = cs[x:x+k]
        fd.write("S\t%s:%s:(A:%s,B:%s)\t%s\n"%(s,x,d[key][0],d[key][1],key))
        lastkmerid[s] = (('%s:%s:(A:%s,B:%s)')%(s,x,d[key][0],d[key][1]))

# Start writing output GFA encoded de-Bruijn Graph
def write_GFA2(G,cs,k,d,fd): 
    global args, g,listofkmers, listoflinks,lastkmerid    
    fd.write("H\tVN:Z:1.0\n")
    for i,x in enumerate(cs):                           # Get the one contig and a number id for the contig
        get_kmers(x,d,k,i,fd)                     # Function to get the fragments of organism A and B if included
    for i,x in enumerate(cs):                           # Get the one contig and a number id for the contig
        get_links(x,d,k,i,fd) 
    # Free up some memory 
    del cs
    for i in G:                                                     # Get the links to the gfa
        for j,o in G[i][0]:
            
            fd.write("L\t%s\t+\t%s\t%s\t%dM\n"%(lastkmerid[i],lastkmerid[j],o,k-1))
        for j,o in G[i][1]:

            fd.write("L\t%s\t-\t%s\t%s\t%dM\n"%(lastkmerid[i],lastkmerid[j],o,k-1))
    # Free up some memory
    del lastkmerid


#########################
#                       #
#  Main Method          #
#                       #
#########################

def main():
    global args
    d = merge_dicts(build(args.A, args.k, args.c), build(args.B, args.k, args.c))
    with open(args.output,"w") as file:
        G,cs = all_contigs(d,args.k)
        write_GFA2(G,cs,args.k,d,file)


parser = argparse.ArgumentParser(description="Creates a GFA file with one or two organisms given a kmer size")
parser.add_argument("-k", default=21, type=int, required=True, help="The kmer size (default = 21)")
parser.add_argument("-A", nargs='+', required=True, help="Organism_A_files")
parser.add_argument("-B", nargs='+', required=True, help="Organism_B_files")
parser.add_argument("-c",  default=0, type=int, required=False, help="Minimum coverage of kmer (default = 0)")
parser.add_argument("-output", default="output.gfa", required=False,help="Output GFA file name")
args = parser.parse_args()
listofkmers = []
listoflinks = []
lastkmerid = {}

# To add more organisms add this parser.add_argument("-B", nargs='+', required=True, help="Organism_B_files")
# change the name and do another call to build and do multiple merge_dicts calls

if __name__ == "__main__":
    main()