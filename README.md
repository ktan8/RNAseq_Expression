RNAseq_Expression
=================

This is a sample code for calculating expression of annotated human genes for an RNA-seq library based on the RPKM methodogy.

Please note some percularities in the code and for situtations where the expression estimates may be a little bit inaccurate. 
Nonetheless, these generally account for only a small proportion of the total reads and the expression estimate for the RPKM 
measurement should be largely inaccurate. In extreme cases, the differences are expected to be atmost in a 10-20% range.

The following code does not account for 

## 1) Multiple splice variants

The Code is unable to extract and distinguish the expression for each and every single splice variants (In fact, most algorithms that claims to be able to do so are not capable of doing it well in every single circumstance more so due to the limitations with RNA-seq datasets rather than the algorithm itself).

This code works by by counting the number of reads that localizes to the annotated locations and then normalizing it to the library size by the RPKM method to give the expression of the annotated transcript.


## 2) Genes with high intronic expressions

In some cases where there are RNA molecules expressed independent and within the intronic regions, the expression of the gene of interest may be a little bit biased. However, such events are not known to occur for most of the genes that we have looked at.



