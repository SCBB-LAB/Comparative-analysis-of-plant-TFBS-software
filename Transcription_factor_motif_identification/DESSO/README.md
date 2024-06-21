# DESSO (DEep Sequence and Shape mOtif) 

DESSO is a deep learning-based framework that can be used to accurately identify both sequence and shape regulatory motifs from the human genome.

## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

Users have the flexibility to choose how they install the necessary packages. However, for efficient package management, we recommend using Anaconda. Once Anaconda is installed, creating and utilizing a virtual environment within Anaconda is a wise option. You can activate a virtual environment with `conda activate` and proceed to install the required packages. If you wish to exit the virtual environment, simply type `conda deactivate`.


#### 1.3 Prerequisites and Dependencies

- Tensorflow 1.1.0 [[Install]](https://www.tensorflow.org/install/)
- CUDA 8.0.44
- Python >= 3.6
- Biopython 1.7.0
- Scikit-learn

To extract the source code for DESSO, execute the following commands:

```
unzip DESSO.zip
```

**Note**

If you want to train your DESSO model on human ChIP-seq data then

- Download [GRCh37.p13.genome.fa](https://bmblx.bmi.osumc.edu/downloadFiles/DESSO/GRCh37.p13.genome.fa.zip) and [encode_101_background](https://bmblx.bmi.osumc.edu/downloadFiles/DESSO/encode_101_background.zip), then unzip them and put them into `data/` directory.
- `data/encode_101`, `data/encode_1001`, and `data/TfbsUniform_hg19_ENCODE` only contain wgEncodeEH002288-related data as an example, owing to the file size limit. To access the source code and whole datasets (totally about 5.9GB) without additional manipulation, just click on [code+whole data](https://bmblx.bmi.osumc.edu/downloadFiles/DESSO/DESSO-master-whole.zip).

## 2. Data information

#### 2.1 Data processing and Model Training Based on Convolutional Neural Network (CNN)
To construct the input datasets for model training kindly follow the instructions here [DESSO](https://github.com/SCBB-LAB/Comparative-analysis-of-plant-TFBS-software/tree/main/DESSO).


#### 3.2 Find motifs on test datset
```
cd code/
python3 predict.py --start_index 0 --end_index 1 --peak_flank 100 --network CNN --feature_format Seq --start_cutoff 0.01 --end_cutoff 1 --step_cutoff 0.03
```
Or simply run the bash script as follows:

```
sh 1.sh
```
**NOTE**: in the bash script `1.sh` kindnly change the TF's name as per your input TF data.

Arguments | Description
----------|----------------------------------------------------------
--start_cutoff | Start of the motif cutoff interval (default is 0.01)
--end_cutoff | End of the motif cutoff interval (default is 1)
--step_cutoff | Increament of the cutoff (default is 0.03)

`--feature_format Seq` indicates that sequence motifs will be predicted. To identify shape motifs, use `--feature_format DNAShape` instead.

**Final Output**
For `--feature_format Seq`, the predicted sequence motifs in different `motif` directory are in `output/encode_101/gc_match/ABF2/Seq/CNN/0`. <br/>
For `--feature_format DNAShape`, four kinds of shape motifs would be predicted.

**Note that** 
if you want to train DESSO with your own negative dataset (non-TFBS sites), then open file `train.py` and open the comment on line numbers 71 and 85 and comment out line 70 and 84.
- Your own dataset in FASTA format with negative dataset must contain four columns with tab separated as shown below.

```
FoldID	Event	seq label
A	peaks	GCGCAAGGCCCATAATATTTTTAGTTATTAAAAAAATTAGCAGACGTAGGGTTGACTTAAAAAAGACTCTTATTACATTAGTCGACAAGTAAAAAACACGTGGCATATATTGTGCGTTCGTAGAGACTGTAATAAAGACGGAGAGATTCTTCTAGAGTCAGTTCTTCTTCTTCATCCTCTTCTTCCCCCCAAATCCTCTCT	1
A	peaks	AACTTTAATTAGTAAAATAGATTTGGCTAAACAAATAAAAAAAACTTTTAGGCTAAAAATTGGATTTGACGTATGAGTAATTGGGGATGAGGGGGACACGTGTCAGAAAATGGGAATGGTATCTTTTGGGGAAAGCATGTAAGTGTGTAATAATGGTCCCCTTCTCTCTCCCATAACCCTACCAAAAATACTTTTCTTTGT	1
A	peak	TGTAAATAAATTGTGTAGCTAATTTGATCTATACAACTATTATTTTTATTAAATATCTATATTTAATCTTATTGTATAAACTTTTTGTTTTACAGCCGACAATTTTTTTTTTTTTTAATATAAAAACATCAGGTTTTGATGAGTGATCTGTTAACAGGGAACGGTCCTACAAAAAGGAACATAGTATACTCTTGATTTTAT	0
A	peak	GAATAGTACGAAAGTAGAGGTGAAACCTTTTTATAATGAAGAGGAAACATTAATTAGCAAGAACCTACATCACATATATTATATATAAGTTCAAACTGCTAAAGATAAAAGTGATTTAATATATACTTGCATTTTTCATTATTAGCAGTCTATCACATGATTCTTTAAGAATAGGTTTGGCTTAGCTAAATTTTTTTTTGG	0
```

Save the `train.py` file and run above mentioned commands for training and prediction with your own negative dataset.


## Citation

If you use DESSO in your research, please cite the following paper:</br>
Jinyu Yang, Anjun Ma, Adam D. Hoppe, Cankun Wang, Yang Li, Chi Zhang, Yan Wang, Bingqiang Liu, and Qin Ma,
"[Prediction of regulatory motifs from human Chip-sequencing data using a deep learning framework](https://academic.oup.com/nar/article/47/15/7809/5542889)",<br/>
Nucleic Acids Research 47, no. 15 (2019): 7809-7824.
