## 1. PeakRanger   
##### Description:  
PeakRanger is a multi-purporse software suite for analyzing next-generation sequencing (NGS) data. 
It contains the following tools:
1. `nr`: a noise ratio estimator useful for QC statistics. Estimates signal to noise ratio which is an indicator for ChIP 
enrichment.
2. `lc`: library complexity calculator useful for QC statistics. Calculates the ratio of unique reads over total reads. 
Only accepts bam files.
3. `ranger`: ChIP-Seq peak caller. Ranger servers better as a narrow-peak caller. It behaves in a conservative but 
sensitive way compared to similar algorithms. It is able to identify enriched genomic regions while at the same time 
discover summits within these regions. 
Ranger supports HTML-based annotation reports.
4. `bcp`: ChIP-Seq peak caller. Tuned for the discovery of broad peaks. BCP supports HTML-based annotation reports.  
5. `ccat`: ChIP-Seq peak caller. Tuned for the discovery of broad peaks. CCAT supports HTML-based annotation reports.
  
Peakranger is installed on [Biowulf.](https://hpc.nih.gov/apps/peakranger.html) 

##### Loading PeakRanger on Biowulf:  

    module load peakranger

##### Running NR (Noise Ratio Estimator):  

    peakranger nr \
    --format bam \
    --data {expt1.bam} \
    --control {control.bam} \
    --output bcp_results
    
##### Running LC (Library Complexity Calculator):  

    peakranger lc \
    --format bam \
    {*.bam} \
    --output bcp_results  

##### Running Ranger (Narrow Peak Caller):  

    peakranger ranger \
    --format bam \
    --report \
    --plot_region 10000 \
    --data {expt1.bam} \
    --control {control.bam} \
    --output bcp_results
    -t 4
 
##### Running BCP (Broad Peak Caller):  

    peakranger bcp \
    --format bam \
    --report \
    --plot_region 10000 \
    --data {expt1.bam} \
    --control {control.bam} \
    --output bcp_results
    -t 4

##### Running CCAT (Broad Peak Caller):  

    peakranger ccat \
    --format bam \
    --report \
    --plot_region 10000 \
    --data {expt1.bam} \
    --control {control.bam} \
    --output bcp_results
    -t 4
<hr>    
