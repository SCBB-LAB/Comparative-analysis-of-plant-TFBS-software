# DeepSTF
## Introduction
DeepSTF, a unique architecture integrating CNN, improved transformer encoder structure, and Bi-LSTM to combine sequence and shape to predict TFBSs

The architecture of the model and the calibration phase steps are explained in this **Figure** from the paper:

<p align="center">
<img src="deepstf.png">
</p>
<p align="center"><b/>Figure: The model workflow</b></p>

## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

Users have their own choice of how to install required packages. But to efficiently manage the installation packages, Anaconda is recommended. After installing Annocoda, it would also be an good option to use virtual environment in annocoda. `conda activate` can be used to activate a virtual environment, and then install required packages. If users want to exit the virtual environment, simply type `conda deactivate`. 

#### 1.2 Install the package and other requirements

**Software list**
- python >=3.6
- pytorch
- numpy 
- pandas
- sklearn
- scipy 
- matplotlib

To extract the source code for DeepSTF, execute the following commands:

```
unzip DeepSTF.zip
```
## 2. Data information

#### 2.1 Data processing
In this part, we will first introduce the **data information** used in this model, then describe the training **data formats**, and finally introduce how to create a data set that meets the model requirements.

Please refer to the example input files **example/ABF2/Sequence/Train_seq.csv** & **example/ABF2/Sequence/Test_seq.csv** in the `example/` directory for train and test files, respectively. The shape files are in the **example/ABF2/Shape/**. To generate the shape files use [DNAshapeR](https://www.bioconductor.org/packages/release/bioc/html/DNAshapeR.html). If you intend to train **DeepSTF** with your own data, make sure to format your data in the same way.

**Note:** Both input files should be in the "csv" format.

## 3. Model Training on integrated CNN, transformer encoder structure, and Bi-LSTM

#### 3.1 Training and evaluation of the DeepSTF model

**Input:** The output files shold be in sequence and shape format for the model training as provided in the `example/ABF2/` directory.

**Usage:**
In the parent directory, execute the following command for model training and the calculation of evaluation metrics on the test dataset:

```
python3.8 train.py ABF2

```
**Output:** 

**Final result:** 

All the trained models from *k*-fold cross-validation, labeled as `params0_bestmodel_(fold_number)fold.hdf5`, and prediction scores (named `score_(fold_number)fold.txt`) from *k*-fold cross-validation are stored in the `output/` directory.

The performance metrics for the test dataset are saved in the `result.txt` file, also located in the `output/` directory.


## Citation

If you use WSCNNLSTM in your research, please cite the following paper:</br>
"[Modeling in-vivo protein-DNA binding by combining multiple-instance learning with a hybrid deep neural network](https://www.nature.com/articles/s41598-019-44966-x)",<br/>
Scientific Reports 9, Article number: 8484 (2019).
