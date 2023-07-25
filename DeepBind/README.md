# DeepBind
## Introduction
This repository contains a PyTorch implementation of a DeepBind model that came out in 2015. The model described in the paper ["Predicting the sequence specificities of DNA and RNA-binding proteins by deep learning"](https://www.nature.com/articles/nbt.3300) is implemented within Pytorch was taken from ["DeepBind-with-PyTorch"](https://github.com/MedChaabane/DeepBind-with-PyTorch). The detailed explanation of the architecture can be found in the [supplementary notes of the paper](https://static-content.springer.com/esm/art%3A10.1038%2Fnbt.3300/MediaObjects/41587_2015_BFnbt3300_MOESM51_ESM.pdf). 

The architecture of the model and the calibration phase steps are explained in this **Figure 1** from the paper:

<p align="center">
<img src="deepbind.jpg">
</p>
<p align="center"><b>Figure: The model workflow</b></p>

## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

Users have their own choice of how to install required packages. But to efficiently manage the installation packages, Anaconda is recommended. After installing Annocoda, it would also be an good option to use virtual environment in annocoda. `conda activate` can be used to activate a virtual environment, and then install required packages. If users want to exit the virtual environment, simply type `conda deactivate`. 

#### 1.2 Install the package and other requirements

Run command to install pytorch

```
python3 -m pip install --pre torch torchvision -f https://download.pytorch.org/whl/nightly/cu111/torch_nightly.html -U
```
Download and extract the source code for DeepBind and move to parent directory, type following commands:

```
unzip DeepBind.zip
```
#### 1.3 Software Requirements

***software list***
- python >=3.6
- pytorch
- numpy 
- pandas
- sklearn
- scipy 
- matplotlib

## 2. Data information

#### 2.1 Data processing

In this part, we will first introduce the **data information** used in this model, then introduce the training **data formats**, and finally introduce how to create a data set that meets to build the model requirements.

We have provided example data format compatible with DeepBind input data (DeepBind input data format: See `example/ABF2_pos.txt`. If you are trying to train DeepBind with your own data, please process your data into the same format as given in above example input data.

## 3. Model Training Based on Convolutional Neural Network (CNN)
#### 3.1 Training and testing 
**Input:** `ABF2_train.txt`, `ABF2_test.txt`. 
All data input files need to be placed in the same folder before training, such as in `example/` directory.

**Usage:**
Run the following command in the parent directory:
```
python3 deepbind.py ABF2
```
**Output:** 

**Final result** 
The trained model and best hyperparameter, `ABF2_Model.pth` and `ABF2_best_hyperpamarameters.pth`, are saved in the `output/` directory, respectively. 
The output file `ABF2_result.txt` located at `output/` directory contains the performance metrics of the test dataset.  

**Note:** A warning massage may appear at the end of the execution of the command as shown below. Please ignore this warning.
```/user_path/python3.8/site-packages/numpy/core/fromnumeric.py:3419: RuntimeWarning: Mean of empty slice.
  return _methods._mean(a, axis=axis, dtype=dtype,
/user_path/python3.8/site-packages/numpy/core/_methods.py:188: RuntimeWarning: invalid value encountered in double_scalars
  ret = ret.dtype.type(ret / rcount)
nan```

## Citation

If you use DeepBind in your research, please cite the following paper:</br>
"[Predicting the sequence specificities of DNA- and RNA-binding proteins by deep learning](https://www.nature.com/articles/nbt.3300)",<br/>
Nature Biotechnology 33, 831–838 (2015).
