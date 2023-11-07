# PlantBind: An attention-based multi-label neural networks for predicting transcription factor binding sites in plants
##   Introduction
The identification of transcription factor binding site(s) (TFBSs) is essential for the analysis of gene regulation. In this context, **PlantBind** is introduced as a method for the integrated prediction and interpretation of TFBS events based on DNA sequences and DNA shape profiles. PlantBind is built upon an attention-based multi-label deep learning framework, enabling it to predict TF binding sites simultaneously. As illustrated in **Figure** , the hybrid model comprises a data processing module, an embedding module, a feature extraction module, and a multi-label output module.

<p align="center">
<img src="flow_chart.png">
</p>  
<p align="center"><b>Figure: The model workflow</b></p>  

The model provides researchers with tools to:

1. Identify the TFBSs of transcription factors.
2. Identify DNA-binding motif of transcription factors

## 1. Environment setup

#### 1.1 Create and activate a new virtual environment

Users have the flexibility to choose how they install the necessary packages. However, for efficient package management, we recommend using Anaconda. Once Anaconda is installed, creating and utilizing a virtual environment within Anaconda is a wise option. You can activate a virtual environment with `conda activate` and proceed to install the required packages. If you wish to exit the virtual environment, simply type `conda deactivate`.

#### 1.2 Install the package and other requirements

Run command to install pytorch

```
python3 -m pip install --pre torch torchvision -f https://download.pytorch.org/whl/nightly/cu111/torch_nightly.html -U
```

Download and extract the source code for PlantBind and move to parent directory, type following commands:


**Software list**
- python      v3.8.12
- pytorch     v1.10.2
- numpy       v1.16.4
- pandas      v1.4.1
- sklearn     v1.0.2
- scipy       v1.5.3
- matplotlib  v3.5.1

To extract the source code for PlantBind, execute the following commands:
```
unzip PlantBind.zip
```

## 2. Data information and processing
In this part, we will first introduce the **data information** used in this model, then introduce the training **data formats**, and finally introduce how to create a data set that meets the model requirements.

All data can be downloaded from the data directory (PlantBind/data)(https://github.com/wenkaiyan-kevin/PlantBind/tree/main/data):
- **Ath-TF-peaks**: the TFBS peak info of 315 Ath TFs, and one neg.bed file
- **Maize-TF-peaks**: the TFBS peak info of 4 maize TFs for trans-species
- **model**: The file that holds the model, which can be loaded to predict new datasets

### 2.1 Construction of the dataset
Next, we will mainly introduce how to create the files mentioned above.
- [Data processing methods](docs/make-datasets.md)
  - [DNA Sequence Data](docs/make-datasets.md#DNA_Sequence_Data)
  - [DNA Shape Data](docs/make-datasets.md#DNA_Shape_Data)
  - [Data Label](docs/make-datasets.md#Data_Label)

#### 2.1.1 DNA Sequence Data

If you intend to train PlantBind with your own data, please ensure that your data is processed into the same format as described above. To begin, you need to split your dataset into training, validation, and test sets. For your convenience, we have provided custom Python programs, namely `dna_shape.py` and `generate_sequence_file.py`, in the parent directory to assist with this conversion process.

Begin by pre-processing your positive and negative FASTA files to partition the data into training, testing, and validation datasets.
```
cat your_pos_fasta_file.fa | paste - -  | awk '{print $2}' > output_file_pos.txt  # For positive fasta file
cat your_neg_fasta_file.fa | paste - - | awk '{print $2}' > output_file_neg.txt  # For negative fasta file
```
To split both the positive and negative data into training, testing, and validation datasets, you can achieve this by executing the following customized Python command:

```
python3 train_valid_test.py output_file_pos.txt output_file_neg.txt train_sequence.table valid_sequence.table test_sequence.table
```
#### 2.1.2 DNA Shape Data
Now, you can use the following customized Python program to convert these ***.table*** files to shape files:
```
python3 dna_shape.py train # For train DNA data shape
python3 dna_shape.py test # For test DNA data shape
python3 dna_shape.py valid # For valid DNA data shape
```
**Output:**
The output files `train_DNA_shape.npyt`, `test_DNA_shape.npy`, and `valid_DNA_shape.npy`, label files are generated.
#### 2.1.3 Data Label
Generating the label matrix is divided into the following steps：  
(i)  Merge all TF peak files into one file  
(ii) Use the `bedtools annotate -i merge.bed -files ... -names ... > model-input-labels.txt` generate the label matrix.  
(iii)Generate binary matrix.
Use following customized python program:
```
python3 generate_binary_matrix.py train 
python3 generate_binary_matrix.py test
python3 generate_binary_matrix.py valid
```
**Output:**
The output files `train_label.txt`, `test_label.txt`, and `valid_label.txt`, label files are generated.

## 3. Model Training Based on attention-based multi-label neural networks
#### 3.1. Data formats for model training and testing

For training, the data mainly consists of three files: (1)DNA sequence file; (2)DNA shape file; (3)data label file  
For predicting, the data mainly consists of three files: (1)DNA sequence file; (2)DNA shape file  

- [Data Format Details Introduction](docs/data_format.md)
  - [DNA Sequence File](docs/data_format.md#DNA_Sequence_File)
  - [DNA Shape File](docs/data_format.md#DNA_Shape_File)
  - [Data Label File](docs/data_format.md#Data_Label_File)


#### 3.2 Training of the model 
 
**Input:** `train_sequence.table`,`train_DNA_shape.npy`,`train_label.txt`.  
All data files need to be placed in the same folder before starting training, such as `data_folder` 

```
python3 -m torch.distributed.launch --nproc_per_node=2 src/main.py --inputs example/ --length 101 --OHEM True --focal_loss True --batch_size 1024 --lr 0.01
```
**Output:** 
**Model** The resulting model file (`trained_model_101_seqs.pkl`) will be saved to `output/` directory.

#### 3.3 Testing of the model 
**Input:** `test_sequence.table`,`test_DNA_shape.npy`,`test_label.txt`.   
All data files need to be placed in the same folder before starting testing, such as `data_folder/` 

```
python main.py --inputs example/ --length 101 --mode test
```
**Output:** 

The output files, including `metrics_b.json`, `metrics_m.json`, `Gmean_threshold_b.txt`, `Gmean_threshold_m.txt`, `precision_recall_curve_TFs_315.pdf`, and `roc_prc_curve_TFs_315.pdf` will be saved to the `output/` directory. Within these output files, `metrics_m.json` file contains information about the evaluation metrics for the test dataset calculated using the multi-label trained model per epoch, while `metrics_b.json` contains information about the evaluation metrics for the test dataset calculated using the binary-label trained model.

## 5. Use the pre-trained model to predict new data

**Input:** `sequence.table`,`DNA_shape.npy`  
The analysis is done by using the [src/inference.py] script and modifying the relevant parameters.

**Output:** 

After training the model, result of the test dataset containing accuracy will be saved to the `output/` directory.

## Citation
If you use PlantBind in your research, please cite the following paper:</br>
[PlantBind: an attention-based multi-label neural network for predicting plant transcription factor binding sites](https://academic.oup.com/bib/article/23/6/bbac425/6713513)",<br/>
Briefings in Bioinformatics 23, no. 6 (2022): bbac425.
