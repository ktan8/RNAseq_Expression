RNAseq_Expression
=================

This is a sample code for calculating expression of annotated human genes for an RNA-seq library based on the RPKM methodogy.

Please note some percularities in the code and for situtations where the expression estimates may be a little bit inaccurate. 
Nonetheless, these generally account for only a small proportion of the total reads and the expression estimate for the RPKM 
measurement should be largely inaccurate. In extreme cases, the differences are expected to be atmost in a 10-20% range.

that the following code does not account for 

## 1) Multiple splice variants

The Code is unable to extract and distinguish the expression for each and every single splice variants (In fact, most algorithms that claims to be able to do so are not capable of doing it well in every single circumstance more so due to the limitations with RNA-seq datasets rather than the algorithm itself).


