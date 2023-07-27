# Comparative-analysis-of-plant-TFBS-software: a repository for TFBS identification tools based on different machine and deep learning tools
This repository provides a hub of all the availble tools that have been developed for TFBS identification for human as well as for plant transcription factors. The corrected and debugged source code links are provided to the evalutated deep learning methods so that the user can easily download and perform the task to identify plant or human TFBS. These the source codes that were developed on various algorithms based on **Machine learning** such as **Support Vector Machine (SVM)**, **XGBoost**, and **deep learning** tools that are based on **Convolutional Neural Network (CNN)**, **Recurrent Neural Network (RNN)**, **Transformer-based** and many other hybrid algorithms. 

## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

Users have their own choice of how to install required packages. But to efficiently manage the installation packages, Anaconda is recommended. After installing Annocoda, it would also be an good option to use virtual environment in annocoda. `conda activate` can be used to activate a virtual environment, and then install required packages. If users want to exit the virtual environment, simply type `conda deactivate`. 

#### 1.2 Install the package and other requirements

Run command to install pytorch

```
python3 -m pip install --pre torch torchvision -f https://download.pytorch.org/whl/nightly/cu111/torch_nightly.html -U
```
Download and extract the source code for Comparative-analysis-of-plant-TFBS-software and move to parent directory, type following commands:

```
git clone https://github.com/SCBB-LAB/Comparative-analysis-of-plant-TFBS-software.git
cd Comparative-analysis-of-plant-TFBS-software
```
Here, you will find all the 14 TFBS identification software repository. Therefore, all you need is to enter into any of the directory and follow the instructions as per the `README.md` of each of the specified software.

## Preparations
If we have 
### Example data

We used *Arabidopsis thaliana* DAP-Seq dataset from [Plant Cistrome Database](http://neomorph.salk.edu/dap_web/pages/browse_table_aj.php). The raw data is derived for 265 plant transcription factors can be easily downloaded.
- Download both [265_dap_data](https://github.com/SCBB-LAB/Comparative-analysis-of-plant-TFBS-software/265_dap_data) of TFBS bed files, then unzip them.
- For generating FASTA data to 
- Here we only provide `ABF2` TFBS data as an example for each of the software, owing to the file size limit.

### Troubleshooting

If there exists any problem in software package installation and module import error, please check the `Supplementary Material S2.docx` from our article for the detailed description of debudding.
If lack other packages when you are running the code, please run `python3.8 pip install -m [package NAME]` directly to the Linux terminal.

If you have any questions/issues/bugs, please post them on [GitHub](https://github.com/SCBB-LAB/Comparative-analysis-of-plant-TFBS-software/issues). They would also be helpful to other users. Or simply drop an email at (dshwljyoti@gmail.com) or (sg927357@gmail.com) if you have any questions.
## Citation
If you use our tool in your research please cite the following paper:</br>
"[Comparative-analysis-of-plant-TFBS-software]",<br/>
Briefings in Bioinformatics.

