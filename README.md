IdCheck
=======
In modern biology research, a growing subdiscipline of genomics is seeking to integrate large-scale 
and different feature data sets in order to achieve a more complete comprehension of a biological 
process. Expression quantitative trait loci (eQTLs) integrated genotype data and gene expression 
data to investigate how genomic loci regulate expression levels of mRNAs or proteins. Before the 
analysis begins, however, it is wise to check whether genotype data and gene expression data come 
from the same sample or not. In this application note, we introduce a tool (IdCheck) for genotype 
and gene expression (RNA sequencing, RNA-seq) sample identity checking. IdCheck compares the 
identity of RNA nucleotide and genotype, and then it calculates a match ratio. Based on this ratio, 
we can determine how the samples are paired. The presented tool provides a fast and efficient way 
to evaluate the potential sample mix problems for tens of geno-type array and RNA-seq lane pairs.