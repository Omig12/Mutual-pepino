# Mutual-pepino

In this repository we will attempt to reproduce a research pipeline that employs a heuristic pairwise alignment algorithm to perform differential expression analysis on Sea Cucumber genetic samples, more specifically on regenerative intestinal tissue samples vs. non regenerative ones.

## Files:
+ **dbg_optimized.py:**
	+ Usage: 
		```sh  
		$ python3 dbg_optimized.py -k <k-mer lenght> -A <fastq A> -B <fastq B> -o <outputfilename> 
		```
	
	+ Ouput:
		+ GFA encoded graph with differential expression information anotated
		```txt  
		H	VN:Z:1.0
		S	0:0:(A:1,B:1)	AAGTTGCGCTAGGGTTAAACT
		S	0:1:(A:1,B:1)	AGTTGCGCTAGGGTTAAACTC
		L	0:0:(A:1,B:1)	+	0:1:(A:1,B:1)	+	20M
		```

+ **analysis.py:**
	+ Usage: 
		```sh 
		$ python3 analysis.py -f <filename> [-t <threshold>] [-o <outputfilename>]
		```

	+ Ouput:
		+ GFA encoded link with differential expression more than or equal to threshold, plus differential expression coeffient.
		```txt  
		L       0:0:(A:1,B:1)   +       0:1:(A:1,B:1)   +       20M
		Differential Coefficient: 0.0
		```

## Worked continued on the [DBGDE](https://github.com/Omig12/DBGDE) repository.

----------------------------------------------------------

## Proposed (DONE) Workflow:
 + Build De-bruijn Graph:
   - Used [this](https://pmelsted.wordpress.com/2013/11/23/naive-python-implementation-of-a-de-bruijn-graph/) approach.
 + Counted the proportion/distribution of kmers recovered from each organism.
 + Perform differential analysis on that data

## Collaborators:
 + Dr. Humberto Ortiz-Zuazaga
 + Kevin Legarreta
 + Louis Gil
 + Israel Dilán-Pantojas
