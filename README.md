# Comparative-analysis-of-plant-TFBS-software

This repository serves as a hub for 32 software tools for peak calling methods, transcription factor regions (TFBRs) identification, and transcription factor binding sites (TFBSs) or motif identification software tools for both human and plant transcription factors. This repository also provides easy to run and executable source codes for all these softwares. Users can easily download, execute, and implement these software tools. The tools covered some of the best performing and most recent machine learning (ML) and deep learning (DL) approaches.

<h2>Machine learning-based algorithms</h2>
<ol>
  <li>Support Vector Machine (SVM)</li>
  <li>XGBoost</li>
</ol> 

<h2>Deep learning-based algorithms</h2>
<ol>
  <li>Convolutional Neural Network (CNN)</li>
  <li>Complex CNN-based algorithms, including: (i) ResNet (ii) DenseNet</li>
  <li>Recurrent Neural Network (RNN)</li>
  <li>Transformers</li>
  <li>Hybrid (i) CNN+RNN and (ii) DesnseNet+Transformer</li>
</ol> 

There are many variations of the ML/DL architectures, as demonstrated in **Figure** below:
 
<p align="center">
<img src="Figure.png">
</p>
<p align="center"><b>Figure: The models' architectures</b></p> 

## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

Users have the freedom to choose how they install the required packages. However, for efficient package management, we recommend using Anaconda. After installing Anaconda, it is also advisable to utilize virtual environments within Anaconda. You can activate a virtual environment using the `conda activate` command and then proceed to install the required packages. If users wish to exit the virtual environment, they can simply type `conda deactivate`. 

#### 1.2 Install the package and other requirements

Run the following command to install PyTorch:

```
python3 -m pip install --pre torch torchvision -f https://download.pytorch.org/whl/nightly/cu111/torch_nightly.html -U
```
Download and extract the source code for Comparative-analysis-of-plant-TFBS-software, and then move to the parent directory. Type the following commands:

```
git clone https://github.com/SCBB-LAB/Comparative-analysis-of-plant-TFBS-software.git
cd Comparative-analysis-of-plant-TFBS-software
```

In this repository, you will discover a collection of 17 TFBS (Transcription Factor Binding Sites) identification software packages. Consequently, all you have to do is navigate into any of these directories, and then follow the instructions outlined in the respective `README.md` file provided for each of the specific software packages.

Each software directory contains its own set of guidelines, usage instructions, and other relevant information in the accompanying `README.md` file.

## 2. Data preparations
If you have an example dataset, the next step is to prepare the training and testing datasets in accordance with the software's input requirements.

### 2.1 DAP-seq data for *A. thaliana*

We used the *Arabidopsis thaliana* DAP-Seq dataset from the [Plant Cistrome Database](http://neomorph.salk.edu/dap_web/pages/browse_table_aj.php). The raw data for 265 plant transcription factors can be easily downloaded for both positive and negative datasets ([from this link 265_dap_data](https://github.com/SCBB-LAB/Comparative-analysis-of-plant-TFBS-software/blob/main/265_dap_data.zip)) in the form of bed files. After downloading, unzip the files. To further process these bed files to FASTA files, follow the instructions described below: 

- For generating FASTA data from bedfiles use either `bedtools getfasta` or `seqtk subseq`.
- For example: `bedtools getfasta -fi genome_file -bed example.bed > out.fa`
- For `seqtk subseq genome_file example.bed > out.fa`

**Note:** Do not forget to save the output FASTA files in their respective text files.

### 2.2 Troubleshooting

If there exists any problem with software package installation or module import error, please check the `Supplementary Material S2.docx` from our article for the detailed description of debugging.
If you lack other packages when you are running the code, please run `python3.8 pip install -m [package NAME]` directly to the Linux terminal.

Post any questions, issues or bugs on the [GitHub repository](https://github.com/SCBB-LAB/Comparative-analysis-of-plant-TFBS-software/issues). They would also be helpful to other users. Or simply drop an email at (dshwljyoti@gmail.com) or (sg927357@gmail.com) if you have any questions.

### Citation
Jyoti, Ritu, Sagar Gupta, Ravi Shankar (2023). Comprehensive evaluation of plant transcription factors binding sites discovery tools. bioRxiv 2023.11.07.566153; doi: https://doi.org/10.1101/2023.11.07.566153 
