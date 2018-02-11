# code taken from https://raw.githubusercontent.com/pmelsted/dbg/master/dbg.py
# import multiprocessing as mp

import argparse
import collections
from Bio import Seq, SeqIO, SeqRecord
import gfapy

#  Create reverse complement
def twin(km):
    return Seq.reverse_complement(km)

# Chunk reads into k sized chunks 
def kmers(seq,k):
    yield (seq[i:i+k] for i in range(len(seq)-k+1))
    # for i in range(len(seq)-k+1):
    #     yield seq[i:i+k]

# Move sequence analysis one nucleotide forward  
def fw(km):
    yield (km[1:]+x for x in 'ACGT')
    # for x in 'ACGT':
    #     yield km[1:]+x

# Move sequence analysis one nucleotide backward
def bw(km):
    # yield(x + km[:-1] for x in 'ACGT')
    for x in 'ACGT':
        yield x + km[:-1]

# Count kmers and build dictionary with counts
# Use limit as cut off for very down regulated seq 
def build(fn,k=31,limit=1):
    d = collections.defaultdict(int)

    # For Each filename in input parse fast.q file 
    for f in fn:
        reads = SeqIO.parse(f,'fastq')
        d[(km for km in (kmers(seq,k) for seq in ((sequence.split('N')) for sequence in ((str(read.seq)) for read in reads))))] += 1
        d[(km for km in (kmers(twin(seq),k) for seq in ((sequence.split('N')) for sequence in ((str(read.seq)) for read in reads))))] += 1
        # d[(km for km in (kmers(twin(segment),k)) for segment in (sequence.split('N')) for sequence in (str(read.seq)) for read in reads)] += 1

        # reads = SeqIO.parse(f,'fastq')
        # # Extract reads 
        # for read in reads:
        # # Maybe it's possible to reduce memory overhead here!    
        #     seq_s = str(read.seq)
        #     seq_l = seq_s.split('N')
        #     # 
        #     for seq in seq_l:
        #         # Extract kmers from forward
        #         for km in kmers(seq,k):
        #             d[km] +=1
        #         # Create Reverse
        #         seq = twin(seq)
        #         # Extract kmers form reverse
        #         for km in kmers(seq,k):
        #             d[km] += 1

 # Remove k-mer elements of dict that appear less than limit times
    (d.delete(i) for i in [x for x in d if d[x] <= limit])
    # d1 = [x for x in d if d[x] <= limit]
    # for x in d1:
    #     del d[x]
    print(d)
    return d

# Merge dictionaries to create homunculus 
def merge_dicts(d1,d2):
    return {i: (d1[i], d2[i]) for i in (list(d1.keys()) + list(d2.keys()))}
    # merge = {}
    # for i in d1.keys():
    #     merge[i] = (d1[i], d2[i])
    # for i in d2.keys():
    #     merge[i] = (d1[i], d2[i])
    # return merge

# Stringify
def contig_to_string(c):
    return (c[0] + ''.join(x[-1] for x in c[1:]))


# we will work here today 
def get_contig(d,km):
    c_fw = get_contig_forward(d,km)
    c_bw = get_contig_forward(d,twin(km))
    # return ((contig_to_string(get_contig_forward(d,km)), get_contig_forward(d,km)) if km in fw(get_contig_forward(d,km)[-1]) else ( (contig_to_string([twin(x) for x in get_contig_forward(d,twin(km))[-1:0:-1]] + get_contig_forward(d,km))), ([twin(x) for x in get_contig_forward(d,twin(km))[-1:0:-1]] + get_contig_forward(d,km))) )

    # return ((contig_to_string(c_fw), c_fw) if km in fw(c_fw[-1]) else ( (contig_to_string([twin(x) for x in c_bw[-1:0:-1]] + c_fw)), ([twin(x) for x in c_bw[-1:0:-1]] + c_fw)) )

    if km in fw(c_fw[-1]):
        c = c_fw
    else:
        c = [twin(x) for x in c_bw[-1:0:-1]] + c_fw
    return contig_to_string(c),c
        
# And here
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

    return c_fw

def all_contigs(d,k):
    done = set()
    r = []
    for x in d:
        if x not in done:
            print(x)
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

    return G,r
    
def get_links(cs,d,k,s,fd):
    for x in range(0,len(cs)-k):
        keyA = cs[x:x+k]
        keyB = cs[(x+1):(x+1)+k]
        
        kmerA = (('%s:%s:(A:%s,B:%s)')%(s,x,d[keyA][0],d[keyA][1]))
        kmerB = (('%s:%s:(A:%s,B:%s)')%(s,x+1,d[keyB][0],d[keyB][1]))

        fd.write("L\t%s\t+\t%s\t+\t%sM\n"%(kmerA,kmerB,(k-1)))

def get_kmers(cs,d, k, s,fd):
    global g, listofkmers,listoflinks,lastkmerid

    for x in range(0,len(cs)+1-k):              #  to get all subsegmet, holds the all the kmers
        key = cs[x:x+k]
        fd.write("S\t%s:%s:(A:%s,B:%s)\t%s\n"%(s,x,d[key][0],d[key][1],key))
        lastkmerid[s] = (('%s:%s:(A:%s,B:%s)')%(s,x,d[key][0],d[key][1]))

def write_GFA2(G,cs,k,d,fd): 
    global args, g,listofkmers, listoflinks,lastkmerid    
    fd.write("H\tVN:Z:1.\n")
    for i,x in enumerate(cs):                           # Get the one contig and a number id for the contig
        get_kmers(x,d,k,i,fd)                     # Function to get the fragments of organism A and B if included
    for i,x in enumerate(cs):                           # Get the one contig and a number id for the contig
        get_links(x,d,k,i,fd) 
    del cs
    for i in G:                                                     # Get the links to the gfa
        for j,o in G[i][0]:
            
            fd.write("L\t%s\t+\t%s\t%s\t%dM\n"%(lastkmerid[i],lastkmerid[j],o,k-1))
        for j,o in G[i][1]:

            fd.write("L\t%s\t-\t%s\t%s\t%dM\n"%(lastkmerid[i],lastkmerid[j],o,k-1))

    del lastkmerid

def main():
    global args
    if args.output:                                             # If the output file name is given use it
        filename = args.output
    else:                                                       # else use standard one
        filename = "output.gfa"

    dA = build(args.A, args.k, 0)
    
    dB = build(args.B, args.k, 0)
    d = merge_dicts(dA, dB)
    del dA
    del dB
    with open(filename,"w") as file:
        G,cs = all_contigs(d,args.k)
        write_GFA2(G,cs,args.k,d,file)


parser = argparse.ArgumentParser(description="Creates a GFA file with one or two organisms given a kmer size")
parser.add_argument("-k", type=int, required=True, help="the kmer size")
parser.add_argument("-A", nargs='+', required=True, help="Organism_A_files")
parser.add_argument("-B", nargs='+', required=True, help="Organism_B_files")
parser.add_argument("-output",required=False,help="Output GFA file name")
args = parser.parse_args()
listofkmers = []
listoflinks = []
lastkmerid = {}

# To add more organisms add this parser.add_argument("-B", nargs='+', required=True, help="Organism_B_files")
# change the name and do another call to build and do multiple merge_dicts calls

if __name__ == "__main__":
    main()