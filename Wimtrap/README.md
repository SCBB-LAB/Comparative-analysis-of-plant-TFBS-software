# Wimtrap
## Introduction
Wimtrap: An integrative tools to predict the location of transcription factor binding sites.

<p align="center">
<img src="wimtrap.jpg">
</p>
<p align="center"><b>Figure: The model workflow</b></p>

## 1. Installation

#### 1.1 Install the package and other requirements

To extract the source code for Wimtrap, execute the following commands:

```
unzip Wimtrap.zip
```
Wimtrap is an R package that requires the last version of R (R 4.0.4), BiocManager and remotes to be installed. 

*Important: the installation of Wimtrap might take up to 1h if all the packages on which it depends need also to be installed.*

In R, type the following lines:
```
if(!require("remotes", quietly = TRUE)){  
    install.packages("remotes")
    }
if(!require("BiocManager", quietly = TRUE)){  
    install.packages("BiocManager")
    }
```
  
Then, you can enter:
```
options(repos = BiocManager::repositories())
getOption("repos")
BiocManager::install("RiviereQuentin/Wimtrap",                     
  dependencies = TRUE,                     
  build_vignettes = TRUE,
  force = TRUE)
quit()    
```

If an error occurs, it might be because the version of R, BiocManager and/or remotes is not updated. 

On linux, it might be also necessary to install as a prerequisite some software. This might be achieved by entering the following in the terminal:

```
sudo apt install libcurl4-gnutls-dev icu-devtools libicu-dev libxml2-dev bzip2-doc libbz2-dev liblzma-dev
```

## 2. Data information

#### 2.1 Data processing

In this section, we will begin by introducing the data information used in this model. Next, we will explain the training data formats, and finally, we will guide you on creating a dataset that aligns with the model's requirements.

We have provided an example data format that is compatible with the Wimtrap input data format (please refer to the `example/` directory).

You can find example input file, namely **ABF2.bed**, in the `example/` directory. If you plan to train Wimtrap with your own data, make sure to prepare your data in the same format.

Before initiating the model training process, you'll need to download a species-specific genome. You can do this by navigating to the `example/` directory. For instance, if you want to download the *A. thaliana* genome:


```
cd example/
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
unzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz
cat Arabidopsis_thaliana.TAIR10.dna.toplevel.fa | cut -d" " -f1 | sed "s/>/>chr/g" | awk '/^>/ {printf("\n%s\n",$0);next;}{printf("%s",$0);} END {printf("\n");}' | grep -A1 ">chr[0-9]" >  GCF_000001735.4_TAIR10.1_genomic-1.fna
cd ../
```

## 3. Model Training Based on XGBoost

#### 3.1 Build and apply a TF-specific model (for example: ABF2)

Predictions can be made by considering chromatin state features related to different conditions. For ***Arabidopsis***, these conditions include whole seedlings, seedling roots, non-hair parts of seedling roots, flowers in stages 4-5, seed coats, heat-shocked seedlings, dark-grown seedlings, dark-grown seedlings exposed to 30 minutes or 3 hours of light, and dark-grown seedlings exposed to a long-day cycle. For ***Solanum***, it encompasses immature and ripening fruits.

To build and apply a TF-specific model, this R package defines specific functions. To predict the binding sites of **ABF2** in Arabidopsis, follow these steps:

```
R
library(Wimtrap)

#The file paths to the genomic data, encoded in BED or GTF/GFF files, are input through the `genomic_data` argument.
#Each file is named according to the feature that it allows to define.
#Remark: the chromosomes are named, in the chrom field of the BED files, according to their number. 
#This number might be preceded by the prefix 'chr' (case-insensitive). For chromosome 1,  'chr1', 'CHR1', 'Chr1' 
#or '1' are accepted.
#Remark: the regions described by a file are all assigned to a score of '1' if the score field is empty (cf. CNS).
# As for the genomic intervals that are not included in a file, they are all assigned to a null score.


args <- commandArgs(trailingOnly = TRUE)
imported_genomic_data.seedlings <- importGenomicData(organism = "Arabidopsis thaliana",
                                                      genomic_data = c(
                                                       DHS = "example/DHS_athal_seedlings_normal.bed",
                                                       DGF = "example/DGF_athal_seedlings_7_days.bed",
                                                       CNS = "example/CNS_athal.bed"
                                                  ))


#The motif representing is encoded in a file in raw pfm format. Other formats are also allowed:
#meme, jaspar, transfac, homer and cis-bp.

#You must specify the name of the transcription factor (here ABF2) as it appears in the file giving
#the motif throught the `TFnames` argument.

#The genome sequence of the considered organism(s) might be automatically downloaded if you provide
#their names through the `organism` argument.
#If you provide the genome sequence from a FASTA file, make sure that the chromosomes are named according to their number. 
#This number might be preceded by the prefix 'chr' (case-insensitive). For chromosome 1,  'chr1', 'CHR1', 'Chr1' 
#or '1' are accepted.


ABF2data.seedlings <- getTFBSdata(pfm = "example/PFMs_athal.pfm",
                               TFnames = "ABF2",
                               organism = "Arabidopsis thaliana",genome_sequence = "example/GCF_000001735.4_TAIR10.1_genomic-1.fna",
                               imported_genomic_data = imported_genomic_data.seedlings)

# Name the `ChIPpeaks` argument according to the training transcription factor(s)

ABF2model <- buildTFBSmodel(ABF2data.seedlings, 
                             ChIPpeaks = c(ABF2 = "example/ABF2.bed"),
                             
                             model_assessment = TRUE)                                
```
**Description of the functions used to build and evaluate a model**
1. importGenomicData() to import the genomic data for the interested species (for example: genomic data for seedlings of Arabidopsis).
2. getTFBSdata() to build the dataset of potential binding sites of specific TF (for example: ABF2) for Arabidopsis thaliana.
3. buildTFBSmodel to get the predictive classifier and, optionally, evaluate the model.
	
**Output:** 

After the model evaluation, three type of output will be generated for test dataset:
1. An annotaion file `annotations.tsv` and model file `model.RData` are saved to parent directory. You may save your output by putting your data by creating TF specific directory.
2. Plot of ROC curve will be saved in the current running directory in `Rplots.pdf` file.  
3. Plot the feature importance, in terms of gain will be saved in the current running directory in `Rplots.pdf` file.
4. Print on the screen the confusion matrix obtained using a prediction score threshold of 0.5, that gives the FP, FN, TP and TN associated to the classification of the potential binding sites of the balanced ‘validation’ dataset by the model.
	
## Citation

If you use Wimtrap in your research, please cite the following paper:<br/>
"[Exploiting Genomic Features to Improve the Prediction of Transcription Factor-Binding Sites in Plants](https://academic.oup.com/pcp/article/63/10/1457/6633738?login=true)",<br/>
Plant and Cell Physiology 63, no. 10 (2022): 1457–1473.
