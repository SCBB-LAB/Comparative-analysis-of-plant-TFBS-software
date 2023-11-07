# TSPTFBS
## Introduction
The architecture of the model and the calibration phase steps are explained in this **Figure** from the paper:

<p align="center">
<img src="TFPTFBS.jpg">
</p>
<p align="center"><b>Figure: The model workflow</b></p>


## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

Users have the flexibility to choose how they install the necessary packages. However, for efficient package management, we recommend using Anaconda. Once Anaconda is installed, creating and utilizing a virtual environment within Anaconda is a wise option. You can activate a virtual environment with `conda activate` and proceed to install the required packages. If you wish to exit the virtual environment, simply type `conda deactivate`.

#### 1.2 Software Requirements
**Software list**

The program requires:
  * python >=3.6
  * tensorflow 2.0.0
  * keras 2.3.1
  * pandas
  * numpy
  * scikit-learn
  * TAIR10 reference genome
  * the [bedtools](https://bedtools.readthedocs.io/en/latest/) software

To extract the source code for TSPTFBS, execute the following commands:
```
unzip TSPTFBS.zip
```

## 2. Data information

#### 2.1 Data processing

In this part, we will first introduce the **data information** used in this model, then introduce the training **data formats**, and finally introduce how to create a data set that meets the model requirements.

We have included an example data format that is compatible with TSPTFBS's input data format (refer to `example/ABF2_pos_train.fa`).

Please review the example input files **ABF2_pos_train.fa** & **ABF2_neg_train.fa** located in the `example/` directory. If you intend to train TSPTFBS with your own data, ensure that your data is prepared in the same format as the provided examples.

## 3. Model Training Based on Convolutional Neural Network (CNN)

#### 3.1 Training TSPTFBS on plant TF datasets
**Input:** `ABF2_pos_train.fa`,`ABF2_neg_train.fa`. 
All data files need to be placed in the same folder before training.

**Note that** both the input files should be in the **FASTA** format.

- **Usage:**
Run the following command in the parent directory:

```
python3 Train.py ABF2
```
**Output:** 

**Final model:** 

The final six trained models, each with a different filter length ranging from 11 to 22, will be saved in the output location `output/ABF2/model/` as `ABF2_pos_train-model-filter_length.hdf5` text files. Additionally, six separate text files, named `ABF2_pos_train-result-filter_length.txt`, based on their respective filter lengths, will be saved in the `output/ABF2/result/` directory. These text files contain the performance metrics for the test dataset.

#### 3.2 Prediction on test dataset 
  
- **Input File Format**

In the `example/` directory, we provide `ABF2_label.txt` and `ABF2_test.fa` files for predicting the test sequences using the pre-trained model. To obtain performance metrics for the test dataset with different filter lengths using the pre-trained model, execute the following command:

```
python3 Predict.py ABF2
```
**Output:**

The prediction results, including accuracy and other metrics (such as the confusion matrix) for the test dataset using various filter_length settings, will also be saved to `ABF2_result.txt`, which can be found in the `output/` directory.

## Citation

If you use TSPTFBS in your research, please cite the following paper:</br>
"[TSPTFBS: a Docker image for trans-species prediction of transcription factor binding sites in plants](https://academic.oup.com/bioinformatics/article/37/2/260/6069568)",<br/>
Bioinformatics 37, no. 2 (2021): 260–262.
