# KEGRU
## Introduction
KEGRU, a model, to identify TF binding sites by combining Bidirectional Gated Recurrent Unit (GRU) network with *k*-mer embedding.

The architecture of the model and the calibration phase steps are explained in this **Figure** from the paper:

<p align="center">
<img src="kegru.jpg">
</p>
<p align="center"><b>Figure: The model workflow</b></p>

## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

Users have the flexibility to choose how they install the necessary packages. However, for efficient package management, we recommend using Anaconda. Once Anaconda is installed, creating and utilizing a virtual environment within Anaconda is a wise option. You can activate a virtual environment with `conda activate` and proceed to install the required packages. If you wish to exit the virtual environment, simply type `conda deactivate`.

#### 1.2 Install the package and other requirements

Install pytorch using following command:

```
python3 -m pip install --pre torch torchvision -f https://download.pytorch.org/whl/nightly/cu111/torch_nightly.html -U
```

**software list**
- python >=3.6
- pytorch
- numpy 
- pandas
- sklearn
- scipy 
- matplotlib

To extract the source code for KEGRU, execute the following commands:

```
unzip KEGRU.zip
```

## 2. Data information

#### 2.1 Data processing
In this part, we will first introduce the **data information** used in this model, then describe the training **data formats**, and finally introduce how to create a data set that meets the model requirements.


We have included an example data format that is compatible with KEGRU's input data format (refer to `example/ABF2_pos.txt`). If you plan to train KEGRU with your own data, ensure that your data is prepared in the same format as described above. It's important to note that the sequences are in *k*-mer format, so you'll need to convert your FASTA format sequences into this format. To facilitate this conversion, we offer a custom Python program named `seq2kmer.py` in the parent directory.

```
python3 seq2kmer.py pos_file.fa 1 > output_file_pos.txt # for positive dataset
python3 seq2kmer.py neg_file.fa 0 > output_file_neg.txt # for negative dataset

```
To split both positive and negative data into training and testing datasets, run the provided customized Python script:

```
python3 train_test.py output_file_pos.tsv output_file_neg.tsv example/ABF2_pos.txt example/ABF2_neg.txt
```
For training dataset: output file `ABF2_pos.txt`
For testing dataset: output file `ABF2_neg.txt`

**Note:** Do not forget to place both files in same directory.

## 3. Model Training Based on Recurrent Neural Network (RNN)

#### 3.1 Training of the model
**Input:** `ABF2_pos.txt`,`ABF2_neg.txt`. 

All data input files should be located in the same folder before training, such as in the `example/` directory. If you intend to train KEGRU with your own data, please ensure that your data is processed into the same format as the provided example.


**Usage:**
Run the following command in the parent directory:

``` 
python3 kegru_train.py -n <file_name> -p <input_file_directory>

Options:

     -g <0-1>     set which gpu (default: 0)
     -n <str>     first word about file name (default: init) (FASTA format)
     -k <int>     set kmer length  (default: 5)
     -s <int>     set stride when slicing k-mers (default: 2)
     -b <float>   set size of one batch (default: 200)
     -i <str>     set initialize vector (default: True)
     -t <str>     set embedding vectors trainable (default: True)
     -l <float>   set the number of layers (default: 1)
     -B <str>     set stride when slicing k-mers (default: True)
     -u <float>   set the number of RNN unit (default: 50)
     -r <str>     result out file
     -O <str>     set hyper parm (default: Adam)
     -N <str>     set hyper parm (default: moname)
     -U <int>     set hyper parm (default: 50)

python3 kegru_train.py -n ABF2 -p example     
```
**Output:** 

**Final result** 

The best trained model and training history files, `ABF2_bestmodel_5_withlstm.hdf5` and `ABF2_training_history.txt`, are saved to the `output/` directory. 

The file named `ABF2_result.txt`, is saved in the `output/` directory, contains the performance metrics for the test dataset.

## Citation

If you use KEGRU in your research, please cite the following paper:</br>
"[Recurrent Neural Network for Predicting Transcription Factor Binding Sites](https://www.nature.com/articles/s41598-018-33321-1)",
Scientific Reports 8, Article number: 15270 (2018).</br>
