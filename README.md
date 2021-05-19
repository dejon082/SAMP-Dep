# SAMP-Dep 
Hackel lab, University of Minnesota
Contact information: dejon082@umn.edu

These Python3 codes analyze deep sequencing read counts to infer antimicrobial peptide activity via quantifying depletion.

Sequence and read count files are labeled by generation, replicate, and inducer concentration. These files must be in the same directory as the analysis codes described below to run properly.

For the first generation, there are two codes, one which calculates slope for individual DNA sequences and one which calculates slope for peptide sequences by combining read counts of synonymous mutants.

For the second generation, there is one code which calculates slope for individual DNA sequences.

For both generations, the output csv file contains the AA or DNA sequence, 1-tailed p-value relative to parental oncocin, average slope (mM^-1 min^-1), average y-intercept (min^-1), slope standard deviation (mM^-1 min^-1), and y-intercept standard deviation (min^-1). 

Output files are included and labeled by generation and DNA versus amino acid sequence based analysis.
