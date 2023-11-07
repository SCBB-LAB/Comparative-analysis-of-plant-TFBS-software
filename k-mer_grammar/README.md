# *k*-mer grammar

## Introduction
The architecture of the model and the calibration phase steps are explained in this **Figure** from the paper:

<p align="center">
<img src="kmer_grammar.jpg">
</p>
<p align="center"><b>Figure: The model workflow</b></p>

## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

Users have the flexibility to choose how they install the necessary packages. However, for efficient package management, we recommend using Anaconda. Once Anaconda is installed, creating and utilizing a virtual environment within Anaconda is a wise option. You can activate a virtual environment with `conda activate` and proceed to install the required packages. If you wish to exit the virtual environment, simply type `conda deactivate`.

#### 1.2 Software Requirements

***software list***
- python >=3.6
- sqlalchemy
- numpy 
- pandas
- sklearn
- scipy 
- matplotlib

To extract the source code for k-mer grammar, execute the following commands:

```
unzip k-mer_grammar.zip
```

## 2. Data information

#### 2.1 Data processing

In this section, we will introduce the data information used in this model, explain the training data formats, and guide you on creating a dataset that aligns with the model's requirements.

We have provided an example data format compatible with *k*-mer grammar input data format (See `example/ABF2_train.txt`).

Please refer to the example input files `ABF2_train.txt` and `ABF2_test.txt` in the `example/` directory. If you intend to train a *k*-mer grammar with your own data, ensure that your data is formatted in the same way.

Each input file, whether for testing or training, consists of two columns separated by a comma. The first column contains `"dna_string"` characters, and the second column contains `"bound"` integers, where 1 indicates a bound state, and 0 indicates a non-bound state.


## 3. Model Training  
#### 3.1 Train and test a "bag-of-*k*-mers" model
To learn how to train a "bag-of-k-mers" model, type the following command:

```
python3 kgrammar_bag-of-k-mer_training_testing.py -help
```
The following command provides information about the algorithm and its usage.

**Usage:**
``` 
python3 kgrammar_bag-of-k-mer_training_testing.py [kmersize: integer] [mode filtered: 'True', or 'False' (i.e., mode full)] [dataset_name: string]
```
**Example:**
``` 
python3 kgrammar_bag-of-k-mer_training_testing.py 8 False ABF2
```

**Output**

**Final result:** 
The above example will train a model with k = 8 without filtering *k*-mers by complexity, reading the file under the `example/` directory. The resulting model file, `kgrammar_bag-of-k-mers_LR_mode_full_ABF2_8_1688192717.pkl`, will be saved to the `output/` directory, along with a database that contains k-mer weights, `kgrammar_bag-of-k-mers_weights_mode_full_ABF2_8_1688192717.db`.

After training the model, the results of the test dataset, including accuracy and other performance metrics (including the confusion matrix), will be saved to `ABF2_grammar_result.txt`. ROC and PRC curve plots, along with the log file, are also located in the `output/` directory.


#### 3.2 Train and test a "vector-*k*-mers" model

To learn how to train a "vector-k-mers" model, type:

```
python3 kgrammar_vector-k-mer_training_testing.py -help
```
To get:

**Usage:**
```
python3 kgrammar_vector-k-mers_training_testing.py [kmersize: integer] [windowsize: integer] [kmer_parsing: 'True' for kmer, 'False' for newtokens] [dataset_name: string]
```
**Example:**
```
python3 kgrammar_vector-k-mer_training_testing.py 8 5 False ABF2
```
**Output**

**Final result:** 

The example above will train a model with k = 8, a window size of 5, and reading the file from the `example/` directory. The vectors for the positive and control sequences will be saved to the `output/` directory, along with the `ABF2_vector_result.txt` file containing performance metrics. After training the model, the script will proceed to test and save ROC and PRC curves, as well as the log file in the same `output/` directory.

## Citation

If you use *k*-mer grammar in your research, please cite the following paper:</br>
"[A *k*-mer grammar analysis to uncover maize regulatory architecture](https://www.nature.com/articles/nbt.3300)",
BMC Plant Biology 19, Article number: 103 (2019).<br/>
