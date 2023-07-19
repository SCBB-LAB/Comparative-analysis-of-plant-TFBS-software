# LS-GKM: A new gkm-SVM software for large-scale datasets

gkm-SVM, a sequence-based method for predicting regulatory DNA elements,
is a useful tool for studying gene regulatory mechanisms.
In continuous efforts to improve the method, new software, `LS-GKM`,
is introduced.  It offers much better scalability and provides further
advanced gapped *k*-mer based kernel functions.  As a result, LS-GKM
achieves considerably higher accuracy than the original gkm-SVM.

gkmSVM-R Tutorial notes

## 1. Installation

After downloading and extracting the source codes and move to parent directory, type:

```

git clone https://github.com/Dongwon-Lee/lsgkm.git
cd lsgkm

```

- Then move to `src/` directory and **make all dependencies** to compile and build source code files.

```   

cd src
make 
make install
```

`make install` will simply copy these executables to the `../bin` directory
    
If successful, you should be able to find the following executables in the current `src/` directory:

    gkmtrain
    gkmpredict
    gkmexplain

Now move back to parent directory
```
cd lsgkm
```

## 2. Tutorial

We introduce the users to the basic workflow of `LS-GKM`.  Please refer to help messages 
for more detailed information of each program.  You can access to it by running the programs 
without any argument/parameter.
  

### 2.1 Training of LS-GKM Based on Support Vector Machine (SVM)
Make sure to make directories to keep your input data and output results before model training.
**Input data directory**
```
mkdir example/
```
**Output data directory**
```
mkdir output/
```
- **Training** 
**Input:** `ABF2_pos_train.fa`,`ABF2_neg_train.fa`, `ABF2_pos_test.fa`,`ABF2_neg_test.fa`. 

All data input files need to be placed in the same folder before training, such as in [example](https://github.com/SCBB-LAB/lsgkm/blob/master/example/). If you are trying to train lsgkm with your own data, please process your data into the same format as it.

Train a SVM classifier using `gkmtrain`. It takes three arguments; 
positive sequence file, negative sequence file, and prefix of output.


    **Usage:** gkmtrain [options] <posfile> <negfile> <outprefix>

     train gkm-SVM using libSVM
```
    Arguments:
     posfile: positive sequence file (FASTA format)
     negfile: negative sequence file (FASTA format)
     outprefix: prefix of output file(s) <outprefix>.model.txt or
                <outprefix>.cvpred.txt
```
    **Options:**
```    
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

First try to train a model using simple test files. Type the following command in parent directory:

```

bin/gkmtrain example/ABF2_pos_train.fa example/ABF2_neg_train.fa output/ABF2_gkmtrain
```

It will generate `ABF2_gkmtrain.model.txt` in the `output/` directory, which will then be used for scoring of 
any DNA sequences as described below.

You can also perform cross-validation (CV) analysis with `-x <N>` option. For example,
the following command will perform 5-fold CV. 

```
bin/gkmtrain -x 5 example/ABF2_pos_train.fa example/ABF2_neg_train.fa output/ABF2_train_gkmtrain
```
The result will be stored in `ABF2_gkmtrain.cvpred.txt` in the `output/` directory.

Please note that it will run SVM training *N* times, which can take time if training 
sets are large.  In this case, you can perform CV analysis on a specific set 
by using `-i <I>` option for parallel runnings. The output will be `<outprefix>.cvpred.<I>.txt`

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
```
    **Options:**
```
     -v <0|1|2|3|4>  set the level of verbosity (default: 2)
                       0 -- error msgs only (ERROR)
                       1 -- warning msgs (WARN)
                       2 -- progress msgs at coarse-grained level (INFO)
                       3 -- progress msgs at fine-grained level (DEBUG)
                       4 -- progress msgs at finer-grained level (TRACE)
    -T <1|4|16>      set the number of threads for parallel calculation, 1, 4, or 16
                     (default: 1)
```
Here, you will try to score the positive and the negative test sequences. Run:

Prediction score for positive sequences. Run following command:

```        
bin/gkmpredict example/ABF2_pos_test.fa output/ABF2_gkmtrain.model.txt output/ABF2_gkmpredict_pos.txt
```
**Output format for postive prediction**

The result will be stored in `output/ABF2_gkmpredict_pos.txt` in the `output/` directory.

Prediction score for negative sequences. Run following command:

```
bin/gkmpredict example/ABF2_neg_test.fa output/ABF2_gkmtrain.model.txt output/ABF2_gkmpredict_neg.txt
```
  
**Output format for negative prediction**

The result will be stored in `output/ABF2_gkmpredict_neg.txt` in the `output/` directory.
  
### 2.3 Evaluating prediction quality 

You may evaluate the model prediction quality as follows:

``` 
python3 scripts/evaluate.py -p output/ABF2
```
    **Options:**
```
    -p 	    prediction file name (prediction score file, output file of gkmpredict for positive and negative file)
```

This will output `output/ABF2_result.txt` all the performance matrics such as accuracy, MCC-values, AUROC, etc.

### 2.4 Generating weight files for deltaSVM

You need to generate all possible non-redundant *k*-mers using the Python script
`scripts/nrkmers.py`.  Then, you score them using `gkmpredict` as described above. 
The output of `lgkmpredict` can be directly used by the deltaSVM script `deltasvm.pl`
available from our deltasvm website.

** Please email Dongwon Lee (dongwon.lee AT childrens DOT harvard DOT edu) if you have any questions. **

## Citation

If you use KEGRU in your research, please cite the following paper:</br>
<br/>
"[Enhanced Regulatory Sequence Prediction Using Gapped k-mer Features](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003711)",<br/>
PLoS Comput Biol 10, e1003711 (2014).
"[LS-GKM: A new gkm-SVM for large-scale Datasets] (https://academic.oup.com/bioinformatics/article/32/14/2196/1742803)",<br/>
Bioinformatics 32, Issue 14, 2196–2198 (2016).
