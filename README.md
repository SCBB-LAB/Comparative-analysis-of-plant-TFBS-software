# Comparative-analysis-of-plant-TFBS-software
This repository provides corrected and debugged source code links to the evalutated deep learning methods from the manuscript on plant TFBS dataset: Evaluation of machine learning and deep learning-based algorithms for discovery of transcription factor binding sites.

The repositories contain the source codes that were developed on various algorithms based on **Machine learning** such as **Support Vector Machine (SVM)**, **XGBoost**, **Convolutional Neural Network (CNN)**, **Recurrent Neural Network (RNN)**, **Transformer-based** and many other hybrid algorithms. 

We tried to evaluate all these software on plant DAP-Seq TFBS dataset. 

## Environment

We evaluated all these software by activating conda environment using following commands:
```
conda activate
python: 3.8.10
pytorch: 1.12.0
torch-geometric: 2.0.1
NVIDIA Driver Version: 525.125.06
CUDA Version: 12.0
GPU: NVIDIA RTX A5000-24Gb
```

## Preparations

### Example data

We used *Arabidopsis thaliana* DAP-Seq dataset from [Plant Cistrome Database](http://neomorph.salk.edu/dap_web/pages/browse_table_aj.php). The raw data is derived for 265 plant transcription factors and are provided with the .
- Download both [265_dap_data](https://github.com/SCBB-LAB/Comparative-analysis-of-plant-TFBS-software/265_dap_data) of TFBS bed or FASTA files, then unzip them.
- Here we only provide `ABF2` TFBS data as an example for each of the software, owing to the file size limit.

### Troubleshooting

If there exists any problem in software package installation and module import error, please check the `Supplementary Material S2.docx` from our article for the detailed description of debudding.
If lack other packages when you are running the code, please run `python3.8 pip install -m [package NAME]` directly to the Linux terminal.

**Please email Jyoti (dshwljyoti@gmail.com) if you have any questions.**
## Citation

If you use Comparative-analysis-of-plant-TFBS-software in your research, please cite the following paper:</br>
<br/>
"[Comparative-analysis-of-plant-TFBS-software]",<br/>
Briefings in Bioinformatics.

