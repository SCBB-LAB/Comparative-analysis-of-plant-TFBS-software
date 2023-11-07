# MAResNet
The architecture of MAResNet are explained in this **Figure 1** from the paper:
<p align="center">
<img src="maresnet.png">
</p>
<p align="center"><b>Figure: The model workflow</b></p>

## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

Users have the flexibility to choose how they install the necessary packages. However, for efficient package management, we recommend using Anaconda. Once Anaconda is installed, creating and utilizing a virtual environment within Anaconda is a wise option. You can activate a virtual environment with `conda activate` and proceed to install the required packages. If you wish to exit the virtual environment, simply type `conda deactivate`. 

#### 1.2 Install the package and other requirements

Run command to install pytorch

```
python3 -m pip install --pre torch torchvision -f https://download.pytorch.org/whl/nightly/cu111/torch_nightly.html -U
```

**Software list**

The required dependencies for MAResNet are in requirements.txt file.

- torchvision==0.9.1      
- pandas==1.2.3
- numpy==1.20.2           
- torch==1.8.1
- scikit_learn==0.24.2

To extract the source code for MAResNet, execute the following commands:
```
unzip maresnet.zip
```

## 2. Data information

#### 2.1 Data processing
In this section, we will first introduce the **data information** used in this model, then explain the training **data formats**, and finally, demonstrate how to create a dataset that adheres to the model's requirements.

We have provided an example data format that is compatible with MAResNet input data format (refer to `example/ABF2/train.data`). If you intend to train MAResNet with your own data, please ensure that your data is processed into the same format and organized into `train.data`, `test.data`, and `valid.data` files.

## 3. Model Training Based on Top-down and bottom-up attentation mechanism and residual network

#### 3.1 Training MAResNet on plant TF datasets
**Input:** `train.data`, `test.data` and `valid.data`.
All data files need to be placed in the same folder before training, such as in the `example/ABF2` directory.

If you want to train this model on your dataset, you should ensure that these three input files are placed in the `example/ABF2/` directory for each TF data. Additionally, make sure that all three input files for different TFs are kept in their respective folders within the `example/` directory.

**Usage:**
To proceed with the next step, please execute the following command in the parent directory:
```
python3 train_on_cell_datasets.py
```

**Output**

**Final result**

The resulting model files, named `maresnet-epoch_number-regular.pth`, will be saved to the `output/checkpoint/ABF2/` directory, alongside regular and best models.
 
After training the model, the result of the test dataset, including accuracy and other metrics (including the confusion matrix) at each epoch, will be saved to a file named `epoch_wise_result.txt` located in the `output/checkpoint/ABF2/` directory. The performance metrics for the best epoch will be output in a file named `test_result.txt` at the same location.

The prediction scores for the test dataset will be saved in `bestiter.pred`, and the log file `df_log2.csv`, which contains performance metrics for the best epoch, will be saved in the `runs/maresnet/ABF2/` directory.

## Citation
If you use MAResNet in your research, please cite the following paper:</br>
<br/>
[MAResNet: predicting transcription factor binding sites by combining multi-scale bottom-up and top-down attention and residual network](https://academic.oup.com/bib/article/23/1/bbab445/6399874)",<br/>
Briefings in Bioinformatics 23, no. 1 (2022).
