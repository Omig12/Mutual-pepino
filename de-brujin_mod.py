# code taken from https://raw.githubusercontent.com/pmelsted/dbg/master/dbg.py

# import multiprocessing as mp
import collections, sys
from Bio import Seq, SeqIO, SeqRecord

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

# count kmers and build dictionary with counts
# use limit as cut off of very down regulated seq 
def build(f,k=31,limit=1):
    d = collections.defaultdict(int)

    # For each filename in input parse fast.q file 
    # for f in fn:
    reads = SeqIO.parse(f,'fastq')
    # Extract reads 
    for read in reads:
        seq_s = str(read.seq)
        seq_l = seq_s.split('N')
        # 
        for seq in seq_l:
            for km in kmers(seq,k):
                d[km] +=1
            seq = twin(seq)
            for km in kmers(seq,k):
                d[km] += 1

 # Remove k-mer elements of dict that
    d1 = [x for x in d if d[x] <= limit]
    for x in d1:
        del d[x]

    return d

def merge_dicts(d1,d2):
    merge = {}
    for i in d1.keys():
        merge[i] = (d1[i], d2[i])
    for i in d2.keys():
        merge[i] = (d1[i], d2[i])
    return merge

# 
def contig_to_string(c):
    return c[0] + ''.join(x[-1] for x in c[1:])


# we will work here today 
def get_contig(d,km):
    c_fw = get_contig_forward(d,km)
    
    c_bw = get_contig_forward(d,twin(km))

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

# Check contigs for kmers and get count
def kmer_count2(cs,d, k, id, file):
    for x in range(0,len(cs)+1-k):
        key = cs[x:x+k]
        # "F 1 * 0 538 1 32 ATGCGCTCGCTCGCTGAGCTGAC A:i:234" 
        file.write("F\t%s\t%d\t%d\t%d\t%d\t%s\tA:i:%s\n"%(id,0,len(cs),x,x+k,key,d[key][0]))
        file.write("F\t%s\t%d\t%d\t%d\t%d\t%s\tB:i:%s\n"%(id,0,len(cs),x,x+k,key,d[key][1]))

# Check contigs for kmers and get count
def kmer_count(cs,d, k, id, file):
    for x in range(0,len(cs)+1-k):
        key = cs[x:x+k]
        # "F 1 * 0 538 1 32 ATGCGCTCGCTCGCTGAGCTGAC A:i:234" 
        file.write("F\t%s\tA:%d\t%d\t%d\t%d\t%s\t*\n"%(id,d[key][0],len(cs),x,x+k,key))
        file.write("F\t%s\tB:%d\t%d\t%d\t%d\t%s\t*\n"%(id,d[key][1],len(cs),x,x+k,key))

# Write to line
def write_GFA2(G,cs,k,d):
    filename = "test.gfa"
    with open(filename, "w+") as file:
        
        file.write("H  VN:Z:2.0\n")
        for i,x in enumerate(cs):
            file.write("S\t%d\t%d\t%s\t\n"%(i,len(x), x ))
            kmer_count(x,d,k,i,file)

        for i in G:
            for j,o in G[i][0]:
                file.write("L\t%d\t+\t%d\t%s\t%dM\n"%(i,j,o,k-1))
            for j,o in G[i][1]:
                file.write("L\t%d\t-\t%d\t%s\t%dM\n"%(i,j,o,k-1))

        for i in d.keys(): 
            file.write("#\t%s\tRC A:%d B:%d\n"%(i,d[i][0],d[i][1])) 

def write_GFA(G,cs,k,d):
    filename = "test.gfa"
    with open(filename, "w+") as file:
        file.write("H  VN:Z:1.0")
        for i,x in enumerate(cs):
            print("S\t%d\t%s\t*\n"%(i,x))
            
        for i in G:
            for j,o in G[i][0]:
                print("L\t%d\t+\t%d\t%s\t%dM\n"%(i,j,o,k-1))
            for j,o in G[i][1]:
                print("L\t%d\t-\t%d\t%s\t%dM\n"%(i,j,o,k-1))

        for i in d.keys(): 
            print("#\t%s\tRC A:%d B:%d"%(i,d[i][0],d[i][1])) 


if __name__ == "__main__":
    if len(sys.argv) < 3: exit("args: <k> <reads_1.fq> ...")
    k = int(sys.argv[1])
    # d = build(sys.argv[2:], k, 1)
    dA = build(sys.argv[2], k, 1)
    dB = build(sys.argv[3], k, 1)
    d = merge_dicts(dA, dB)
    with open('dict_of_ratio_1.txt', "w+") as f:
        f.write(str(d))
    G,cs = all_contigs(d,k)
    write_GFA2(G,cs,k,d)