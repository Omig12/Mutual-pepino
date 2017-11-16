# Mutual-pepino

In this repository we will attempt to reproduce a research pipeline that employs a heuristic pairwise alignment algorithm to perform differential expression analysis on sea cucumber samples, more specifically on regenerative intestinal tissue samples vs. non regenerative ones.

## Workflow:
 + Build De-bruijn Graph
   - Might use [this](https://pmelsted.wordpress.com/2013/11/23/naive-python-implementation-of-a-de-bruijn-graph/) approach.
 + Count the proportion/distribution of kmers recovered from each organism.
 + Perform differential analysis on that data
