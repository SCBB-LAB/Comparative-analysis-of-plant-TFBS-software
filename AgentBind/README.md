# AgentBind
## Introduction

<a href="https://zenodo.org/badge/latestdoi/174050946"><img src="https://zenodo.org/badge/174050946.svg" alt="DOI"></a>

AgentBind is a machine-learning framework for analyzing context regions of binding sites and identifying specific non-coding nucleotides with strong effects on binding activities. This code repository contains code for the classification + visualization experiments with the DanQ and DeepSEA architectures respectively.

## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

All experiments are executed on CentOS Linux 7 (core) with Python (v2.7.5). Prior to your code execution, please make sure you have installed the following tools/libraries.

Our code requires external python libraries including tensorflow v1.9.0 GPU-version, biopython v1.71, numpy v1.15.4, six v1.14.0, scikit-image v0.14.5, and matplotlib. You can install them with the pip package manager:

```pip install numpy six matplotlib biopython sklearn scikit-image tensorflow-gpu==1.9.0```


#### 1.2 Install FIMO from the MEME-suite

Download and install `MEME-suite` software if not installed earlear by using folloing command line:

1)    Download the software from https://meme-suite.org/doc/download.html
2)    Type the following commands
      - tar zxf meme-5.5.6.tar.gz
      - cd meme-5.5.6
      - ./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
      - make
      - make test
      - make install
            
3)    Edit your shell configuration file to add $HOME/meme/bin and $HOME/meme/libexec/meme-5.5.6 to your shell's path. This can often be done by editing the file named .profile to add the following line:
    export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.6:$PATH

## 2. Data information

#### 2.1 Data processing

**Data for experiments with the DanQ architecture**
https://drive.google.com/file/d/12mrLk9Ci7u2tKB8kuqldGXE9ghAzpbUk/view?usp=sharing

**Data for experiments with the DeepSEA architecture**
https://drive.google.com/file/d/1UaaqgFlce9FSaBX2RoIz9pDaXacwQ3lW/view?usp=sharing

If the datasets are not downloaded from above mentioned links, then try to construct your own datasets with following instructions and command lines:

It is up to the user to choose between the `DanQ` or `DeepSea` architecture. Therefore, they should navigate to the directory and follow the instructions.

To preprocess your ChIP-seq dataset, ensure that the identified motifs are in Position Weight Matrix (PWM) format, as shown in example/ABF2/meme.txt. This file will be used to scan the motif in the genome file of the specified species (e.g., Arabidopsis thaliana).

To run the motif scanning tool `FIMO` from the `MEME-suite`, use the following command in your terminal:

```
fimo --oc example/ABF2_fimo example/ABF2/meme.txt example/genomes/tair10/tair10.fa
```

This will generate the `example/ABF2_fimo` directory containing a `fimo.tsv` file and other summary files. Before proceeding, remove the last four lines from the `fimo.tsv` file to avoid an IndexError:

```
python3 data_pre_software.py  --fimo_file example/ABF2_fimo/fimo.tsv --chipseq example/ABF2_narrow.bed --scope 1000 --resultdir example/ABF2 --datadir example --blockcore c --recorddir example/ABF2
```

If you encounter an error regarding the missing `hg19.fa.fai` file, ensure that the genome file for your species is present in the `genomes/` directory (e.g., tair10.fa.fai for A. thaliana). Additionally, modify `data_prep.py` at lines 52, 53, 256, and 257, and `converter.py` at line 55.

After these adjustments, run the following command:

```
python3 data_prep/data_prep.py  --fimo_file example/ABF2_fimo/fimo.tsv --chipseq example/ABF2_narrow.bed --scope 1000 --resultdir example/ABF2 --datadir example --blockcore c --recorddir example/ABF2
```

This command will generate three files of `training`, `test`, and `validation` files within `example/ABF2/` directory, which consists of train, test, and validation datasets for a particular TF dataset. These three files are further used to train, test, and validate the model. No run following command to train the model:

#### 2.2 Model Training Based on Recurrent Neural Network (RNN)

- **Training** 
**Input:** `example/ABF2/training/data.txt`,`example/ABF2/validation/data.txt`. 

All data input files need to be placed in the same folder before training, such as in `example/`. If you are trying to train KEGRU with your own data, please process your data into the same format as it.

**Usage:**
Run the following command in parent directory:


```
python2.7 ../model/train-transfer-learning.py  --data_dir example/ABF2/training --valid_dir example/ABF2/validation --n_train_samples 8974 --n_valid_samples 2238 --train_dir ABF2_out --seq_size 1000 --checkpoint_dir ABF2_out  --batch_size 128 --n_classes 2
```
The output checkpoints and trained models are saved to `ABF2_out/` directory.


## Run ##
**AgentBind.py** is the go-to python script which execute all the experiments.

**Required parameters:**
* --datadir: the directory where you stored the downloaded data.
* --motif: a text file containing the names, motifs, and ChIPseq files of TFs of interest. This text file can be found in the given data at `{your-data-path}/table_matrix/table_core_motifs.txt`.
* --workdir: a directory where to store all the intermediate/oversized files including the well-trained models, one-hot-encoded input sequences, and Grad-CAM annoation scores.
* --resultdir: a directory where to store all the results.

To run AgentBind, you can simply execute:
```
python AgentBind.py 
--motif {your-data-path}/table_matrix/table_core_motifs.txt 
--workdir {your-work-path}
--datadir {your-data-path}
--resultdir {your-result-path}
```

AgentBind reports results of two situations, core motifs present (c) and blocked (b). You can find the correspounding classification results (AUC curves) in: `{your-result-path}/{b or c}/{TF-name}+GM12878/`. And the Grad-CAM annoation scores are available at `{your-work-path}/{TF-name}+GM12878/seqs_one_hot_{b or c}/vis-weights-total/weight.txt`.

The python program "AgentBind.py" takes ~24-48 hours to complete. If you need the Grad-CAM annotation scores only, you can directly download them here (DanQ version only):
* https://drive.google.com/file/d/1HB-_bG1K6rbbtBxh2OQ2ldL5BVp3NlBQ/view?usp=sharing

For questions on usage, please open an issue, submit a pull request, or contact Melissa Gymrek (mgymrek@ucsd.edu) or An Zheng (anz023@eng.ucsd.edu).

## Citation

Published paper: [**click here**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8009085/)

Please cite: \
`Zheng, A., Lamkin, M., Zhao, H. et al. Deep neural networks identify sequence context features predictive of transcription factor binding. Nat Mach Intell (2021).`

