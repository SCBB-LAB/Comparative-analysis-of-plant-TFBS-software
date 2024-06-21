Introduction
============
PeakRanger is a multi-purporse software suite for analyzing next-generation sequencing (NGS) data. 
The suite contains the following tools:

1. `nr` noise rate estimator. Estimates signal to noise ratio which is an indicator for ChIP enrichment
2. `lc` library complexity calculator. Calculates the ratio of unique reads over total reads. Only accepts bam files.
3. `wig` coverage file generator. Generates variable step format wiggle file
4. `wigpe` coverage file generator. Generates bedGraph format wiggle file and supports spliced alignments and thus only supports bam files
5. `ranger` ChIP-Seq peak caller. It is able to identify enriched genomic regions while at the same time discover summits within these regions.
6. `ccat` ChIP-Seq peak caller. Tuned for the discovery of broad peaks
7. `bcp` ChIP-Seq peak caller. Tuned for the discovery of broad peaks. 

`ranger`, `ccat` and `bcp` support HTML-based annotation reports.

`wigpe` can also generate coverage files for bam files containing spliced reads, such as those from RNA-Seq experiments.

If you use PeakRanger in your research, please cite:

Feng X, Grossman R, Stein L: PeakRanger:A cloud-enabled peak caller for ChIP-seq data.BMC Bioinformatics 2011, 12(1):139.(http://www.biomedcentral.com/1471-2105/12/139/)

If you use the `ccat` tool, please also cite:

Xu, H., L. Handoko, et al. (2010).A signal-noise model for significance analysis of ChIP-seq with negative control.Bioinformatics 26(9): 1199-1204.(http://bioinformatics.oxfordjournals.org/content/26/9/1199)

If you use the `bcp` tool, please also cite:

Xing H, Mo Y, Liao W, Zhang MQ (2012) Genome-Wide Localization of Protein-DNA Binding and Histone Modification by a Bayesian Change-Point Method with ChIP-seq Data. PLoS Comput Biol 8(7): e1002613. doi:10.1371/journal.pcbi.1002613.(http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002613)

System Requirement
===================
1. The ranger and ccat tool depends on the R programming environment to generate HTML reports. But they can run without it, even with `--report` and `--gene_annot_file` set.
2. When the number of peaks called by ranger or ccat is huge, it takes a while for user's browser to parse the generated HTML file.
3. The lc tool needs about 1.7G ram per 10 million aligned reads.

How PeakRanger works 
=================

##ranger
 Calling narrow peaks
It turns out that ranger servers better as a narrow-peak caller. It behaves in a conservative but sensitive way compared to similar algorithms.

 The algorithm
ranger uses a staged algorithm to discover enriched regions and the summits within them. In the first step, PeakRanger implements a FDR based adapative thresholding algorithm, which was originally proposed by PeakSeq. PeakRanger uses this thresholder to find regions with enriched reads that exceed expects. After that, PeakRanger searches for summits in these regions. The summit-search algorithm first looks for the location with largest number of reads. It then searchs for sub-summits with the sensitivity, the delta -r, specified by the user. Smaller -r will generate more summits.The coverage profiles are smoothed and padded before calling summits. The smoothing grade varies with -b. Higher smoothing bandwidth results less false summits at the cost of degraded summit accuracy .To measure the significance of the enriched regions, PeakRanger uses binormial distribution to model the relative enrichment of sample over control. A p value is generated as a result. Users can thus select highly significant peaks by using a smaller -p.

 Reads extending
ranger extends reads before calling peaks. The default reads extension length is 200. However, users can change this by -l if the datasets come with a different fragment size. The extension length will change the reads coverages generated from the raw reads as it will change the heights of peaks.

##wigpe and wig
To help visualizing the results, wigpe and wig generates reads coverage files in the wig format. These files can then be loaded into browsers to evaluate the authenticity of called peaks. Since smaller wiggle files take less time and memory to load, --split can be set to generate one small wig file per chromosome.

##ccat
 Calling broad peaks
Calling broad peaks remain unsolved for the ChIP-Seq community. It seems the CCAT algorithm is one of those that is designed for this problem, especially for calling histone modification marks.

 The algorithm
For details of the algorithm, please refer to the original manuscript of CCAT:

Xu, H., L. Handoko, et al. (2010).A signal-noise model for significance analysis of ChIP-seq with negative control.Bioinformatics 26(9): 1199-1204.(http://bioinformatics.oxfordjournals.org/content/26/9/1199)

##bcp
 Calling broad peaks
bcp also serves as a broad peak caller. In many situations we perfer bcp over ccat while there are certain scenarios ccat outperforms. A drawback of the current implementation is that it does not support summit calling which is supported by `ranger` and `ccat`.

 The algorithm
For details of the algorithm, please refer to the original manuscript of bcp:

Xing H, Mo Y, Liao W, Zhang MQ (2012) Genome-Wide Localization of Protein-DNA Binding and Histone Modification by a Bayesian Change-Point Method with ChIP-seq Data. PLoS Comput Biol 8(7): e1002613. doi:10.1371/journal.pcbi.1002613.(http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002613)

##nr
nr is a module of the original CCAT algorithm that estimates the similarity of data and control. It indicates roughly how data departs from control

##lc
lc measures the percentage of unique reads. The result measures how diversified the reads are in the dataset. The idea is from:

Chen, Yiwen, Nicolas Negre, Qunhua Li, Joanna O. Mieczkowska, Matthew Slattery, Tao Liu, Yong Zhang, et al. 2012. Systematic evaluation of factors influencing ChIP-seq fidelity. Nature Methods 9(6): 609-614.(http://www.nature.com/nmeth/journal/v9/n6/full/nmeth.1985.html)

Compiling PeakRanger from source codes
======================================

Compiling in Ubuntu
------------------

Required libraries before compiling:

1. The Boost library v1.56 or newer

2. Pthread

3. g++

4. zlib

Install any missing items via `apt-get install` .
Once all the libraries are installed, go to the root path of the unzipped package and type:

 `make`

This will generate `bin/peakranger`. Compilation in other Linux distributions is similar.

We suggest that you modify the `STATIC` variable in the makefile:

  `STATIC = -static`

This will produce a stand-alone binary file that is easier to use.

Compiling in Mac OSX
------------------

Required libraries before compiling:

1. Xcode developer tool kit from Apple

2. The Boost library v1.56 or newer

3. zlib

4. Pthread

The Xcode kit can be installed using the OSX installation disk. If you dont have the installation disk, you can also get it for free from Apple Developer. The tool kit installs essential command line tools such as make and C++ compilers. The Boost library can be installed by following the instructions on its website. If you do not have root access,
modify the BOOST_PATH variable in the make file:

  `BOOST_PATH = -I/path/to/your/boost/header -L/path/to/your/boost/library`

Once all the libraries are installed, go to the root path of the unzipped package and type:

  `make` 

If the compilation failed, double check the BOOST_PATH variable is correctly set.
The resulting binaries require dynamic boost library files, to make sure `peakranger` can find
these files, type :

  `export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/path/to/your/boost/library`

please change the `path` accordingly.

Compiling in Windows
------------------
Not supported but should be possible.

Synopsis
===================
  peakranger lc sample.bam 
  
  peakranger nr --format bam sample.bam control.bam
  
  peakranger wig --format bam sample.bam sample.bam_coverage
  
  peakranger wig --format bed sample.bed sample.bed_coverage
  
  peakranger wigpe sample.bam sample.bam_coverage 
  
  peakranger wigpe sample.bam sample.bam_coverage_splitted -s
  
  peakranger wigpe sample.bam sample.bam_coverage_splitted_by_strand -sx
  
  peakranger wigpe sample.bam sample.bam_coverage_gzipped -z
  
  peakranger wigpe sample.bam sample.bam_coverage_splitted_by_strand_gzip -sxz
  
  peakranger wigpe sample.bam sample.bam_coverage_read_extend_to_200 -l 200
  
  peakranger ranger --format bam sample.bam control.bam ranger_result 
  
  peakranger ranger --format bam sample.bam control.bam ranger_result_threaded_faster -t 3
  
  peakranger ccat --format bam sample.bam contro.bam ccat_result
  
  peakranger ccat --format bam sample.bam contro.bam ccat_result_with_HTML_report --report --gene_annot_file hg19refGene.txt 
                  
  peakranger ccat --format bam sample.bam contro.bam ccat_result_with_HTML_report_5kb_region --report --gene_annot_file hg19refGene.txt --plot_region 10000
  
  peakranger bcp --format bam sample.bam contro.bam bcp_result 
  
  peakranger bcp --format bam sample.bam contro.bam bcp_result_with_HTML_report_5kb_region --report --gene_annot_file hg19refGene.txt --plot_region 10000

Command line options
===================
## nr
  input

  `-d,--data`
  
data file.

  `-c,--control`	
  
control file. 
  
  `--format`
  	
the format of the data file, can be one of : bowtie, sam, bam and bed.

  Qualities

  `-l,--ext_length`   

read extension length

  Other 

  `-h,--help`   
  
  show the usage

  `--verbose`  
  
  show progress

  `--version`  
  
  output the version number

## lc
  input

  `-d,--data`
  
data file.

  Other 

  `-h,--help`   
  
  show the usage

  `--verbose`  
  
  show progress

  `--version`  
  
  output the version number

## wig
  input

  `-d,--data`
  
data file.

  `--format`
  	
the format of the data file, can be one of : bowtie, sam, bam and bed.

  Output

  `-o,--output`  
  
  the output location

  `-s,--split`
        
 generate one wig file per chromosome
  
  `-z,--gzip`
        
 compress the output
  
  `-x,--strand`
        
 generate one wig file per strand
  
  Qualities

  `-l,--ext_length`   

read extension length

  Other 

  `-h,--help`   
  
  show the usage

  `--verbose`  
  
  show progress

  `--version`  
  
  output the version number

## wigpe
  input

  `-d,--data`
  
data file.

  Output

  `-o,--output`  
  
  the output location

  `-s,--split`
        
 generate one wig file per chromosome
  
  `-z,--gzip`
        
 compress the output
  
  `-x,--strand`
        
 generate one wig file per strand
  
  Qualities

  `-l,--ext_length`   

read extension length

  Other 

  `-h,--help`   
  
  show the usage

  `--verbose`  
  
  show progress

  `--version`  
  
  output the version number

## bcp
  input

  `-d,--data`
  
data file.

  `-c,--control`	
  
control file. 
  
  `--format`
  	
the format of the data file, can be one of : bowtie, sam, bam and bed.

  Output

  `-o,--output`  
  
  the output location

  `--report`
        
 generate html reports
  
  `--plot_region`
        
 the length of the snapshort regions in the HTML report. It also controls the search span for nearby genes.
  
  `--gene_annot_file`
        
 the gene annotation file
  
  Qualities

  `-p,--pval`   
    
p value cut off

  `--win_size`   
    
Sliding window size

  `-l,--ext_length`   

read extension length

  Other 

  `-h,--help`   
  
  show the usage

  `--verbose`  
  
  show progress

  `--version`  
  
  output the version number

## ranger
  input

  `-d,--data`
  
data file.

  `-c,--control`	
  
control file. 
  
  `--format`
  	
the format of the data file, can be one of : bowtie, sam, bam and bed.

  Output

  `-o,--output`  
  
  the output location

  `--report`
        
 generate html reports
  
  `--plot_region`
        
 the length of the snapshort regions in the HTML report. It also controls the search span for nearby genes.
  
  `--gene_annot_file`
        
 the gene annotation file
  
  Qualities

  `-p,--pval`   
    
p value cut off

  `-q,--FDR`   
    
FDR cut off

  `-l,--ext_length`   

read extension length

  `-r,--delta`   

sensitivity of the summit detector

  `-b,--bandwidth`    
  
smoothing bandwidth.

  `--pad`    

pad read coverage to avoid false positive summits

  Running modes

  `-t`  
      
number of threads.(default: 1)
  
  Other 

  `-h,--help`   
  
  show the usage

  `--verbose`  
  
  show progress

  `--version`  
  
  output the version number

## ccat
  input

  `-d,--data`
  
data file.

  `-c,--control`	
  
control file. 
  
  `--format`
  	
the format of the data file, can be one of : bowtie, sam, bam and bed.

  Output

  `-o,--output`  
  
  the output location

  `--report`
        
 generate html reports
  
  `--plot_region`
        
 the length of the snapshort regions in the HTML report. It also controls the search span for nearby genes.
  
  `--gene_annot_file`
        
 the gene annotation file
  
  Qualities

  `-q,--FDR`   
    
FDR cut off

  `--win_size`   

 sliding window size

  `--win_step`   

window moving step

  `--min_count`    
  
 minimum window reads count

  `--min_score`    

 minimum window reads fold change

  `-l,--ext_length`   

read extension length

  Other 

  `-h,--help`   
  
  show the usage

  `--verbose`  
  
  show progress

  `--version`  
  
  output the version number

Output files
===============
  lc 
lc does not generate any files

  nr
nr does not generate any files

  wig
wig generates a single file by default. When -x or -s is specified, it generates multiple files depending on the datasets.

  wigpe
similar to wig

  ranger 
Three files will be geneated:

`_summit.bed`

`_region.bed`

`_details`

The first two bed files can be visualized in IGV. `_summit.bed` file contains
the locations of summits ranked by their FDR. `_regions.bed` file contains the locations
of regions ranked by their FDR. Each summit or region is annotated by the 4th column.

`_details` file contains both summits and regions as well
as the regions's FDR and p values. 

When `--report` is enabled, this file will also contain nearby genes of called peaks.

`--report` enables HTML reporting that generates a folder named using the data file's name. The folder
contains a single `index.html` visualizable in most browsers.

  ccat 
Similar to ranger

  bcp
Similar to ranger but it does not provide the `_summit.bed` file.

PeakRanger: http://ranger.sourceforge.net
modENCODE: http://www.modencode.org
PeakSeq: http://info.gersteinlab.org/PeakSeq
sourceforge: http://ranger.sourceforge.net
ChIP-Seq: http://en.wikipedia.org/wiki/Chip-Sequencing
GNU Scientific Library (GSL): http://www.gnu.org/software/gsl/
Ubuntu: http://www.ubuntu.com
confuse option file parser: http://www.nongnu.org/confuse/
Bowtie: http://bowtie-bio.sourceforge.net
Hadoop: http://hadoop.apache.org/
Boost: http://www.boost.org/
IGV: http://www.broadinstitute.org/igv/
Hadoop Streaming: http://hadoop.apache.org/common/docs/r0.15.2/streaming.html
NGS: http://en.wikipedia.org/wiki/DNA_sequencing
