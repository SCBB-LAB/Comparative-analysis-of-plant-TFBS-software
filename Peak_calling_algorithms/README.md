# Benchmarking ChIP-Seq Peak Callers

## Table of Contents
> 1. **[PeakRanger](#1-peakranger)**
> 2. **[PePr](#2-pepr)**
> 3. **[MUSIC](#3-music)**
> 4. **[Ritornello](#4-ritornello)**
> 5. **[MACS2](#5-macs2)**
> 5. **[LanceOtron](#6-lanceotron)**
<hr>

## 1. PeakRanger   
##### Description:  
PeakRanger is a multi-purporse software suite for analyzing next-generation sequencing (NGS) data. 
It contains the following tools:
1. `ranger`: ChIP-Seq peak caller. Ranger servers better as a narrow-peak caller. It behaves in a conservative but sensitive way compared to similar algorithms. It is able to identify enriched genomic regions while at the same time discover summits within these regions. Ranger supports HTML-based annotation reports.
2. `lc`: library complexity calculator useful for QC statistics. Calculates the ratio of unique reads over total reads. Only accepts bam files.
3. `nr`: a noise ratio estimator useful for QC statistics. Estimates signal to noise ratio which is an indicator for ChIP enrichment.
4. `bcp`: ChIP-Seq peak caller. Tuned for the discovery of broad peaks. BCP supports HTML-based annotation reports.
5. `ccat`: ChIP-Seq peak caller. Tuned for the discovery of broad peaks. CCAT supports HTML-based annotation reports.

##### Peak calling using PeakRanger 

```
./peakranger ranger -d input.bam -c control.bam -o out_file.txt --format bam
```

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

**NOTE:** For adding more options to your command line:

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

<hr>

## 2. PePr
##### Description:
PePr is a ChIP-Seq Peak-calling and Prioritization pipeline that uses a sliding window approach and models read counts across replicates and between groups with a negative binomial distribution. PePr empirically estimates the optimal shift/fragment size and sliding window width, and estimates dispersion from the local genomic area. Regions with less variability across replicates are ranked more favorably than regions with greater variability. Optional post-processing steps are also made available to filter out peaks not exhibiting the expected shift size and/or to narrow the width of peaks. 


##### Running Pepr:
```
PePr -c chip_rep1.bam,chip_rep2.bam -i input_rep1.bam,input_rep2.bam -f bam -s 100 --peaktype=sharp --output-directory=pepr_out
```

**NOTE:** For multiple samples if one sample for input file then we need to provide only one input data file.
    
* --shiftsize Half the fragment size. This should be half of the size that we used for macs
<hr> 

## 3. MUSIC
##### Description: 
MUSIC is a tool for identification of enriched regions at multiple scales in the read depth signals from ChIP-Seq experiments. 

##### Download and Installation

You can download MUSIC C++ code <a href="https://github.com/gersteinlab/MUSIC/archive/master.zip">here</a>. There are no dependencies for building MUSIC. After download, type:

unzip MUSIC.zip<br>
cd MUSIC<br>
make clean<br>
make
</font></i>
</div><br>
to build MUSIC. The executable is located under directory <font face="courier">bin/</font>. It may be useful to install <a href="http://samtools.sourceforge.net/">samtools</a> for processing BAM files.

To get help on which options are available, use:
<div style="padding:8px;background-color:#ddd;line-height:1.4;">
<font face="courier">
MUSIC -help
</font>
</div>
##### Running Music:
    mkdir chip; mkdir input
    samtools view chip.bam | MUSIC -preprocess SAM stdin chip/ 
    samtools view input.bam | MUSIC -preprocess SAM stdin input/
    samtools view /directory/to/chip.bam | MUSIC -preprocess SAM stdin chip/ 
    samtools view /directory/to/input.bam | MUSIC -preprocess SAM stdin input/
    mkdir chip/sorted;mkdir chip/dedup;mkdir input/sorted;mkdir input/dedup
    MUSIC -sort_reads chip chip/sorted 
    MUSIC -sort_reads input input/sorted 
    MUSIC -remove_duplicates chip/sorted 2 chip/dedup 
    MUSIC -remove_duplicates input/sorted 2 input/dedup
      
    MUSIC -get_multiscale_broad_ERs \
    -chip chip/dedup \
    -control input/dedup \
    -mapp Mappability_36bp \
    -l_mapp 36 \
    -begin_l 1000 \
    -end_l 16000 \
    -step 1.5
   * This code tells MUSIC to identify the enriched regions starting from 1kb smoothing window length upto 16kb with 
   multiplicative factor of 1.5 using the default parameters for the remaining parameters. The ERs for each scale are dumped.
   
* You will have to generate appropriate mappability files (https://github.com/gersteinlab/MUSIC#multi-mappability-profile-generation) Our read lengths are 100 or 125, so generate for those. Also, make sure you have bowtie2 indices for the genomes. Look at https://github.com/gersteinlab/MUSIC#running-music-with-default-parameters-and-automatic-selection-of-l_p-parameter- to automate parameter selection.*
<hr>


## 4. Ritornello 
Ritornello is a ChIP-seq peak calling algorithm based on signal processing that can accurately call binding events without the need to do a pair total DNA input or IgG control sample.  It has been tested for use with narrow binding events such as transcription factor ChIP-seq.

[Ritornello Preprint](http://biorxiv.org/content/early/2015/12/11/034090)

##### Download precompiled binaries

Currently compiled for ubuntu x64. Mac and Windows versions coming soon.

[Download](https://github.com/kstant0725/Ritornello/releases)

##### Compiling on ubuntu:

install dependencies:

-samtools

`sudo apt-get install libbam-dev`

-FFTW

`sudo apt-get install libfftw3-dev`

-boost

`sudo apt-get install libboost-dev`

Checkout the source using:

`git clone https://github.com/KlugerLab/Ritornello.git`

Compile:

`cd Ritornello`

`make`

The executable is then made at `bin/Ritornello`.
You can move it where ever you like and/or add it to your path, or simply run it from its bin location.

#### Creating a sorted bam file:

This tutorial assumes the user starts with sequenced ChIP-seq reads in the fastq format, `MyFile.fastq` for singled end or `MyFile_1.fastq` and `MyFile_2.fastq` for paired-end.

The first step in preparing a sorted bam file, the input for Ritornello, is to map the ChIP-seq reads to the reference genome. This can be done using any comparative genome alignment tool.  For this tutorial we will use Bowtie, which can be downloaded [here](http://bowtie-bio.sourceforge.net/index.shtml)

Also install samtools, which can be downloaded [here](http://samtools.sourceforge.net/) or on Ubuntu installed using:

`sudo apt-get install samtools`

Map the reads to the genome using the following command.  In this example hg19 (human genome version 19) is the prefix for the bowtie index (reference genome), so make sure to download the correct index (also found on the bowtie website for most organisms) and specify the prefix that is appropriate for your organism.  `MyFile.fastq` for singled end or `MyFile_1.fastq` and `MyFile_2.fastq` for paired-end are user provided fastq files containing the ChIP-seq reads from the sequencer. For single end run:

`bowtie -S -n 2 -k 1 -m 1 -X 1000 -I 0 --best --strata hg19 MyFile.fastq | samtools view -bS - > MyBamFile.bam`

or if you have paired end data, run:

`bowtie -S -n 2 -k 1 -m 1 -X 1000 -I 0 --best --strata -1 MyFile_1.fastq -2 MyFile_2.fastq hg19 | samtools view -bS - > MyBamFile.bam`
        
This will create a bam (alignment) file in the same directory.  See the bowtie documentation for a detailed explanation of each option.

To create a sorted bam file, call the samtools sort command on your bam file as follows:

`samtools sort MyBamFile.bam MySortedBamFile`

Indexing the bam file is useful so that it can be used with other tools such as the IGV genome browser, but not required to run ritornello.  To index the bam file, run the following:

`samtools index MySortedBamFile.bam`

This will create an index `MySortedBamFile.bam.bai` file in the same directory

#### Using Ritornello:
-basic usage

`./Ritornello -f MySortedBamFile.bam -o file_name`

Where `MySortedBamFile.bam` is an index/sorted bam file that can be obtained by first mapping the fastq files using an aligner (such as bowtie) and then sorting and indexing using samtools.


The full script with details on how to create all required files and run the script is provided
[here](https://github.com/KlugerLab/Ritornello/blob/master/Scripts/AnalyzeRitornelloOutput.R)
 
#### Ritornello options:
`--help`	print the help message

`--version`	print Ritornello's current version

`-f <ChIP.bam>`	ChIP.bam is a ChIP-seq sorted bam file to call peaks on.  If you're bamfile is not sorted, please use the samtools sort utility.  Additionally, running the samtools index utility may be required by some visualization tools (IGV etc.) but is not required by Ritornello.

`-o <OutputPrefix>`	Specifies a prefix to use when reporting output.  This can be a file path.
Ex. `-o /home/MyUser/MyOutputPrefix` would report called peaks to `/home/MyUser/MyOutputPrefix-peakSummary.narrowPeak`

`-p <int>`	maximum number of threads to use during execution

`-q <decimal>`	The -log10(q-value) used to threshold reported peaks.  Ex.  Specifying `-q 2` (the default) will report all peaks that are more significant than q-value=0.01.  Set `-q 0` when calling peaks for input to the Irreproducible Discovery Rate software.

`-s <decimal>`	The signal value used to threshold reported peaks.  The signal value is the effect size for the reported peak and has units of read count.  It is the sum of beta1 for the positive and negative strands around a reported peak, and can best be interpreted as the maximum likelihood number of reads due to binding at the reported peak.  Ex. specifying `-s 40` (the default) will report all peaks with signal value greater than 40.  Set `-s 0` when calling peaks for input to the Irreproducible Discovery Rate software.

`-n <decimal>`	The minimum read and matched filter threshold.  This specifies the minimum number of reads per window of size twice the maximum fragment length centered around the peak required to perform a likelihood ratio test.  Additionally the matched filter (which is also normalized to units of read counts is thresholded by this number. This threshold is mostly used to control compute time by limiting the number of tests performed.  `-n 20` is the default.  Setting it too high (larger than effect size) may cause lower expressed peaks to be missed.  Setting it lower generally increases runtime and memory usage. Set `-n 10` when calling peaks for input to the Irreproducible Discovery Rate software.

`--OCE`	Specifying the OCE option tells Ritornello to include an additional term in the likelihood ratio test to more strictly control for Open Chromatin Effects.  These are areas of high coverage which are generally uniform and also present in sonicated input DNA.  Ritornello's background coverage term can usually control for most open chromatin effects, however, when coverage is extremely high, it can call spurious peaks at the bounderies of these regions where coverage changes abruptly.  `--OCE` is useful to avoid spurious results in highly sequenced small genomes (yeast), but may cause a loss of sensitivity and not recommended for mouse, human, etc.

`--no-artifact-handling`	Specifying the `--no-artifact-handling` option tells Ritornello not to try to detect and remove read length artifacts.  You can add this option if your reads are fairly long and you suspect there wont be any issues with artifacts related to mismapping of reads.

#### Ritornello advanced options:

`--debug-folder`	Specify a folder to print out debug files (must end with a "/").  This is mainly used for development purposes.

`--filter-file`	Specify a filter shape to use for this run.  Filter shapes are printed to fir.txt in the debug folder when it is specified.  It is simply a vector of numbers giving the negative strand (or reverse positive strand) filter shape.  --FLD-file must be set to use this option.  Mainly used for development purposes

`--FLD-file`		Specify a fragment length distribution for use with this run.  FLD files are printed to fld.txt in the debug folder when the option is specified.  It is simply a vector giving the distibution of the fragment lengths.  Mainly used for development purposes

<hr>

## 2. MACS2 
##### Description:
MACS empirically models the length of the sequenced ChIP fragments, which tends to be shorter than sonication or library construction size estimates, and uses it to improve the spatial resolution of predicted binding sites. MACS also uses a dynamic Poisson distribution to effectively capture local biases in the genome sequence, allowing for more sensitive and robust prediction. MACS compares favorably to existing ChIP-Seq peak-finding algorithms and can be used for ChIP-Seq with or without control samples.


##### Running MACS2 (Narrow Peak Mode):  

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

## 6. LanceOtron 
##### Description:

**LanceOtron** is a machine learning, genomic data extraction and analysis tool trained for ATAC-seq, ChIP-seq, and DNase-seq peak calling. A freely available and fully-featured webtool version, utilising the graphical user interface [MLV](https://mlv.molbiol.ox.ac.uk) and hosted at the [MRC WIMM Centre of Computational Biology, University of Oxford](https://www.imm.ox.ac.uk/research/units-and-centres/mrc-wimm-centre-for-computational-biology), can be found at [LanceOtron.molbiol.ox.ac.uk](https://lanceotron.molbiol.ox.ac.uk).

## Python Package and Requirements

LanceOtron was built using Python 3.8.3 and TensorFlow 2. The models have been saved such that a TensorFlow 1 setup could be used making only minor amendments to the scripts (see note in modules folder). Additional packages were used for benchmarking LanceOtron - see [requirements.txt](lanceotron/requirements.txt) for specific version numbers used.

LanceOTron is currently tested on Linux and MacOSX, please see the CircleCI implementation for more information if interested. We do not currently support native Windows installations, nor do we have any immediate plans due to the amount of work involved to maintain it. Windows subsystem for Linux should function as expected, though if you have issues using LanceOTron on your specific installation, please [raise an issue](https://github.com/LHentges/LanceOtron/issues/new/choose) and we will help you resolve them. 

**Additional Python Packages for Benchmarking:**

N.B. [bedtools](https://github.com/arq5x/bedtools2) needs to be installed to use *pybedtools*.
> * pandas
> * matplotlib
> * pybedtools
> * seaborn

## Command Line Installation

There are three ways to install LanceOTron. The first and second methods are recommended for general users while the second is recommended for developers interested in extending the source code.

### Method 1: Installing from Pypi 

LanceOTron is [hosted on Pypi](https://pypi.org/project/lanceotron/) and can be easily installed using `pip`. 

```sh
pip install lanceotron
```

### Method 2: Conda installation

Conda installation is available via [sgriva's](https://anaconda.org/sgriva) channel.

```sh
conda install -c sgriva lanceotron
```

### Method 3: Local installation

We recommend using a fresh virtual environment with Python 3.7+. 

1. Clone the repository and navigate to the CLI package
2. Install dependencies with pip.
3. Install the package.
4. Run tests to ensure that everything is working.

```{sh}
git clone git@github.com:LHentges/LanceOtron.git # Step 1
cd LanceOTron/lanceotron

pip install -r requirements.txt # Step 2

pip install -e . # Step 3

python -m unittest # Step 4
```

## Usage

Currently there are 3 LanceOtron modules available. By default the modules return a bed file of **all candidate peaks, good, bad and otherwise**, along with their associated scores. For more details regarding candidate peak selection and the deep neural network scoring please see the citation below. 

Detailed usage instructions for each of the three modules along with a tutorial are available in the CLI directory [here](lanceotron/).

### Preparing bigwig files

All LanceOtron modules require a bigwig file to supply the model with coverage data. We recommend directly converting BAM files to bigwigs with [deepTools](https://github.com/deeptools/deepTools/tree/develop) using the following command:

> `bamCoverage --bam filename.bam.sorted -o filename.bw --extendReads -bs 1 --normalizeUsing RPKM`

The options used in this command are important, as they affect the shape of peaks and therefore the neural network's assessment. Extending the reads out to the fragment length represents a more accurate picture of coverage (N.B for paired end sequencing the extension length is automatically determined, single end tracks will require the user to specify the `--extendReads` length), as does using a bin size of 1 (the `--bs` flag). We recommend RPKM normalisation, as this was also used for the training data.

## Citation

Please see our [Bioinformatics article](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btac525/6648462) for further details on this project or if you use it in your own work.

```bibtex
@article{10.1093/bioinformatics/btac525,
    author = {Hentges, Lance D and Sergeant, Martin J and Cole, Christopher B and Downes, Damien J and Hughes, Jim R and Taylor, Stephen},
    title = "{LanceOtron: a deep learning peak caller for genome sequencing experiments}",
    journal = {Bioinformatics},
    year = {2022},
    month = {07},
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btac525},
    url = {https://doi.org/10.1093/bioinformatics/btac525},
    note = {btac525},
    eprint = {https://academic.oup.com/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btac525/45048211/btac525.pdf},
}
```

<hr>
