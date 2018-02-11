# Test Space

In this folder you will find:
+ seq_test.fq
	+ The First sequence of an E. coli FastQ file 
	+ 70 nucleotides long 

+ seq2_test.fq
	+ A piece of the same sequence of an E. coli FastQ file 
	+ 10 nucleotides long 

+ seq_test.gfa
	+ Results of running optimized_dbg.py on seq_test.txt against itself
	+ $ python3 -k 10 -A seq_test.fq -B seq_test.fq -output seq_test.gfa
		* The result was as expected:
			* (70-k)+1 k-mers found
			* All equally expressed
			* Completely Connected graph 


+ seq_test2.gfa
	+ Results of running optimized_dbg.py on seq_test.txt against seq2_test.txt
	+ $ python3 -k 10 -A seq_test.fq -B seq2_test.fq -output seq_test2.gfa
		* The result was as expected: 
		* (70-k)+1 k-mers found
		* Only the first k-mer was shared
		* Completely Connected graph 