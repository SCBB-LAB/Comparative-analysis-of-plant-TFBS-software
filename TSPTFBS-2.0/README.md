# TSPTFBS-2.0
TSPTFBS 2.0 is a [webserver](http://www.hzau-hulab.com/TSPTFBS/) based on deep learning models for transcription factor binding site (TFBS) prediction. It can be used to mine the potential core motifs within a given sequence by the trained 389 TFBS prediction models of three species (Zea mays, Arabidopsis, Oryza sativa) and the three interpretability algorithms Deeplift，*in-silico* tiling-deletion and *in-silico* mutagenesis. [TSPTFBS 2.0](https://github.com/liulifenyf/TSPTFBS-2.0) is freely accessible for all users. 

## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

Users have the flexibility to choose how they install the necessary packages. However, for efficient package management, we recommend using Anaconda. Once Anaconda is installed, creating and utilizing a virtual environment within Anaconda is a wise option. You can activate a virtual environment with `conda activate` and proceed to install the required packages. If you wish to exit the virtual environment, simply type `conda deactivate`.

#### 1.2 Install the package and other requirements

Run following command to install pytorch:

```
python3 -m pip install --pre torch torchvision -f https://download.pytorch.org/whl/nightly/cu111/torch_nightly.html -U
```

The program requires:
  * python==3.7.13
  * tensorflow-gpu==2.0.0 (for model training)
  * tensorfow==1.14.0 (for DeepLIFT, TF-MoDISco and model predicting)
  * deeplift
  * modisco
  * keras==2.3.1.0.2
  * scikit-learn==1.0.2
  * pandas 
  * numpy 
  * the [bedtools](https://bedtools.readthedocs.io/en/latest/) software
To extract the source code for LS-GKM, execute the following commands:

```
unzip TSPTFBS-2.0.zip

```

###  Predicting (389 models were employed to predict the binding intensity of inpuy sequences)
```
cd TSPTFBS-2.0
python predict.py <input fasta file>
```
After running the program, a file named 'results.csv' will be generated in the current folder which records the prediction results of the models.
We here provide a test.fa file for an example: 
```
python predict.py Example/test.fa
```
#### Input File Format
The input file must contain DNA sequences which have a length of 500bp with a FASTA format.
A FASTA file of the example is:
```
>4:175156999-175157499
GAATGTGCGTGCTGTGTTGCAGTCGCGTTAGGGCCAAGTCCTAGCCTTTGTGGTGATTAGATTTAGGGGGTGGTCAAGATTCACATATTTATGTTTCTTAACCCTCTCCTGGACTTGGCGACTCTTTTTTTTACCCCCTCCCGAGACAAGTGCCCGTGCGTTTCTTGTTGAACTCTGAATTTGCTTATTCAACAGAAGTTGATAATGATAATAAAAGAAGAGGCATCCTGTGTAAATCGATGCCTCATTTTCTTACTGCCTGTCAGGCTGTCATGGCATGTCAGCAGCTGGGACGGAGATTTGCATGTAAATGTTGTACAGAATTGCATGATCTATCCTGTGAAGCAGAATCAAAATTCTGCTCGGGTAAGATAATGATAAACAGCATAGATGCTGGCTATATGTGTACGAGTACTTGCTACAAAGTGAACCATGGAGCACTTTCTTTTTGATAATTACCATGGTGCAGGTTGAGATGCGAGAATGTTGTATGCCGAGAC
```
#### Output File Format
The output file will seem like below: the first column represents the names of 389 TFs, the remaining columns (The example has one remaining column because the input file has one enquired DNA sequences) record the probabilities of given DNA sequences to be predicted as a TFBS of one of 389 TFs.
```
TF  4:175156999-175157499	
AT3G10113	5.8268368E-05 
AT3G12130	0.0003848466	
AT3G52440	0.6031477	
...
```
### interpretability
```
cd TSPTFBS-2.0
python interpretability.py <input fasta file> <species> <tf>
```
It should be noted ```<species>``` that one is chosen from 'Zea_mays_models','Arabidopsis_models' and 'Oryza_sativa_models'.
It should be noted ```<tf>``` that one is chosen from the tf names of selected species.
After running the program, a dir about deeplift results will be generated in the current folder.
We here provide a test.fa file and employed one of models of Zea mays for an example: 
```
python interpretability.py Example/test.fa Zea_mays_models ALF2 
```
### TF-MoDISco
```
cd TSPTFBS-2.0
python modisco_test.py <input fasta file> <species> <tf>
```
It should be noted ```<species>``` that one is chosen from 'Zea_mays_models','Arabidopsis_models' and 'Oryza_sativa_models'.
It should be noted``` <tf>``` that one is chosen from the tf names of selected species.
After running the program, a dir about tf-modisco results will be generated in the current folder.
We here provide a test.fa file and employed one of models of Zea mays for an example: 
```
python modisco_test.py Example/test.fa Zea_mays_models ALF2 
```
### Citation
* Huang, G., et al. Densely Connected Convolutional Networks. IEEE Computer Society 2016.
* Shrikumar, A., Greenside, P. and Kundaje, A. Learning Important Features Through Propagating Activation Differences. 2017.
* Shrikumar, A., et al. Technical Note on Transcription Factor Motif Discovery from Importance Scores (TF-MoDISco) version 0.5.6.5. In.; 2018. p. arXiv:1811.00416.
