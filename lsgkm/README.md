# LS-GKM: A new gkm-SVM software for large-scale datasets
## Introduction
gkm-SVM, a sequence-based method for predicting regulatory DNA elements, is a useful tool for studying gene regulatory mechanisms. In continuous efforts to improve the method, new software, `LS-GKM`, is introduced.  It offers much better scalability and provides further advanced gapped *k*-mer based kernel functions.  As a result, LS-GKM achieves considerably higher accuracy than the original gkm-SVM.

## 1. Installation

To extract the source code for LS-GKM, execute the following commands:

```
unzip lsgkm.zip
```
- After navigating to the `src/` directory, you can proceed to build the source code and its dependencies by using the following command:

```    
cd src
make 
```
    
If successful, you should be able to find the following executables in the current `src/` directory:

    gkmtrain
    gkmpredict

Then run in the linux terminal:
```    
make install
```
`make install` will simply copy these executables to the `../bin` directory

Now move back to the parent directory
```
cd ../
```

## 2. Tutorial

To access more detailed information about each program, please refer to the help messages. You can access these help messages by running the programs without providing any arguments or parameters. These help messages typically contain valuable information about the program's options and usage.
  

### 2.1 Training of LS-GKM Based on Support Vector Machine (SVM)

- **Training** 
**Input:** `ABF2_pos_train.fa`,`ABF2_neg_train.fa`, `ABF2_pos_test.fa`,`ABF2_neg_test.fa`. 

Before training, ensure that all data input files are placed in the same folder, such as the `example/` directory. If you intend to train LS-GKM with your own data, you should process your data to match the required format used in the `example/` directory. This ensures that the training process runs smoothly and effectively.

Train a SVM classifier using `gkmtrain`. It takes three arguments; positive sequence file, negative sequence file, and prefix of output.


    **Usage:** gkmtrain [options] <posfile> <negfile> <outprefix>

```
    Arguments:
     posfile: positive sequence file (FASTA format)
     negfile: negative sequence file (FASTA format)
     outprefix: prefix of output file(s) <outprefix>.model.txt or <outprefix>.cvpred.txt

    Options:
 
     -t <0 ~ 5>   set kernel function (default: 4 wgkm)
                  NOTE: RBF kernels (3 and 5) work best with -c 10 -g 2
                    0 -- gapped-kmer
                    1 -- estimated l-mer with full filter
                    2 -- estimated l-mer with truncated filter (gkm)
                    3 -- gkm + RBF (gkmrbf)
                    4 -- gkm + center weighted (wgkm)
                         [weight = max(M, floor(M*exp(-ln(2)*D/H)+1))]
                    5 -- gkm + center weighted + RBF (wgkmrbf)
     -l <int>     set word length, 3<=l<=12 (default: 11)
     -k <int>     set number of informative column, k<=l (default: 7)
     -d <int>     set maximum number of mismatches to consider, d<=4 (default: 3)
     -g <float>   set gamma for RBF kernel. -t 3 or 5 only (default: 1.0)
     -M <int>     set the initial value (M) of the exponential decay function
                  for wgkm-kernels. max=255, -t 4 or 5 only (default: 50)
     -H <float>   set the half-life parameter (H) that is the distance (D) required
                  to fall to half of its initial value in the exponential decay
                  function for wgkm-kernels. -t 4 or 5 only (default: 50)
     -R           if set, reverse-complement is not considered as the same feature
     -c <float>   set the regularization parameter SVM-C (default: 1.0)
     -e <float>   set the precision parameter epsilon (default: 0.001)
     -w <float>   set the parameter SVM-C to w*C for the positive set (default: 1.0)
     -m <float>   set cache memory size in MB (default: 100.0)
                  NOTE: Large cache signifcantly reduces runtime. >4Gb is recommended
     -s           if set, use the shrinking heuristics
     -x <int>     set N-fold cross validation mode (default: no cross validation)
     -i <int>     run i-th cross validation only 1<=i<=ncv (default: all)
     -r <int>     set random seed for shuffling in cross validation mode (default: 1)
     -v <0 ~ 4>   set the level of verbosity (default: 2)
                    0 -- error msgs only (ERROR)
                    1 -- warning msgs (WARN)
                    2 -- progress msgs at coarse-grained level (INFO)
                    3 -- progress msgs at fine-grained level (DEBUG)
                    4 -- progress msgs at finer-grained level (TRACE)
    -T <1|4|16>   set the number of threads for parallel calculation, 1, 4, or 16
                     (default: 1)
```


To begin, start by training a model using simple test files. You can initiate the process by entering the following command in the parent directory:


```
bin/gkmtrain example/ABF2_pos_train.fa example/ABF2_neg_train.fa output/ABF2_gkmtrain
```

After executing the training command, it will generate `ABF2_gkmtrain.model.txt` in the `output/` directory. This model file will be used for scoring of any DNA sequences as described below.

Additionally, you can also perform cross-validation (CV) analysis with the `-x <N>` option. For example, the following command will perform 5-fold CV. 

```
bin/gkmtrain -x 5 example/ABF2_pos_train.fa example/ABF2_neg_train.fa output/ABF2_train_gkmtrain
```
The results will be saved in `ABF2_gkmtrain.cvpred.txt` in the output/ directory.

Please note that SVM training will be run *N* times, which can be time-consuming for large training sets. In such cases, you can perform CV analysis on a specific set by using the `-i <I>` option for parallel runs. The output will be saved as `<outprefix>.cvpred.<I>.txt`.

The format of the cvpred file is as follows:
  
    [sequenceid] [SVM score] [label] [CV-set]
    ...

### 2.2 Scoring DNA sequence using gkm-SVM

You use `gkmpredict` to score any set of sequences.

    **Usage:** gkmpredict [options] <test_seqfile> <model_file> <output_file>

     score test sequences using trained gkm-SVM
```
    Arguments:
     test_seqfile: sequence file for test (FASTA format)
     model_file: output of gkmtrain
     output_file: name of output file

Options:

     -v <0|1|2|3|4>  set the level of verbosity (default: 2)
                       0 -- error msgs only (ERROR)
                       1 -- warning msgs (WARN)
                       2 -- progress msgs at coarse-grained level (INFO)
                       3 -- progress msgs at fine-grained level (DEBUG)
                       4 -- progress msgs at finer-grained level (TRACE)
    -T <1|4|16>      set the number of threads for parallel calculation, 1, 4, or 16
                     (default: 1)
```
Here, you will calculate prediction scores for both the positive and negative test sequences.

- To get the prediction scores the positive sequences, please execute the following command:

```        
bin/gkmpredict example/ABF2_pos_test.fa output/ABF2_gkmtrain.model.txt output/ABF2_gkmpredict_pos.txt
```
**Output format for positive prediction**

The result will be stored in `output/ABF2_gkmpredict_pos.txt` file in the `output/` directory.

- To obtain the prediction scores for the negative sequences, please execute the following command:

```
bin/gkmpredict example/ABF2_neg_test.fa output/ABF2_gkmtrain.model.txt output/ABF2_gkmpredict_neg.txt
```
  
**Output format for negative prediction**

The results will be stored in the `output/ABF2_gkmpredict_neg.txt` file in the `output/` directory.
  
### 2.3 Evaluating prediction quality 

You may evaluate the model prediction quality as follows:

``` 
python3 scripts/evaluate.py -p output/ABF2

Options:

    -p 	    prediction file name (prediction score file, output file of gkmpredict for positive and negative file)
```

This command will output the performance metrics, including accuracy, MCC-values, AUROC, etc., in the `output/ABF2_result.txt` file.


### 2.4 Generating weight files for deltaSVM

You need to generate all possible non-redundant *k*-mers using the Python script `scripts/nrkmers.py`. Then, you score them using `gkmpredict` as described above. The output of `lgkmpredict` can be directly used by the deltaSVM script `deltasvm.pl` available from our deltasvm website.

## Citation

If you use LS-GKM in your research, please cite the following paper:</br>
<br/>
"[Enhanced Regulatory Sequence Prediction Using Gapped k-mer Features](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003711)",
PLoS Comput Biol 10, e1003711 (2014).
"[LS-GKM: A new gkm-SVM for large-scale Datasets] (https://academic.oup.com/bioinformatics/article/32/14/2196/1742803)",
Bioinformatics 32, no. 14 (2016): 2196–2198. </br>
