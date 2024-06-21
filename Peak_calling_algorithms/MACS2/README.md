# MACS2 

## Description:
MACS empirically models the length of the sequenced ChIP fragments, which tends to be shorter than sonication or library construction size estimates, and uses it to improve the spatial resolution of predicted binding sites. MACS also uses a dynamic Poisson distribution to effectively capture local biases in the genome sequence, allowing for more sensitive and robust prediction. MACS compares favorably to existing ChIP-Seq peak-finding algorithms and can be used for ChIP-Seq with or without control samples.


## Running MACS2 (Narrow Peak Mode):  

macs2 callpeak -t TF.bam -c control.bam --format=BAM --name=TF --gsize=genome_size_of_specise --tsize=26

Input file options

-    -t: The IP data file (this is the only REQUIRED parameter for MACS)
-   -c: The control or mock data file
-    -f: format of input file; Default is “AUTO” which will allow MACS to decide the format automatically.
-    -g: mappable genome size which is defined as the genome size which can be sequenced; some precompiled values provided.

Output arguments

-   --outdir: MACS2 will save all output files into speficied folder for this option
-   -n: The prefix string for output files
-   -B/--bdg: store the fragment pileup, control lambda, -log10pvalue and -log10qvalue scores in bedGraph files

<hr>
