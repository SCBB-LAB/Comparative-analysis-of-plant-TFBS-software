# Introduction
Contrast Motif Finder (CMF) is a de novo tool for discovering differentially enriched motifs in two sets of sequences, offering non-discretized PWM estimations corrected for false positives. It generates motifs predominantly found in either experimental or control sequences and provides binding site lists for both sets. For each motif, CMF calculates a likelihood ratio (LR) score for identified sites. To run CMF on your input data, enter into the paraent directory:
```
cd CMF
```

## 1. Data information

#### 1.1 Data processing

In this section, we will begin by introducing the **data information** utilized in method. Next, we will describe the training **data formats**. Finally, we will provide instructions on how to create a dataset that adheres to the model's requirements. For the input dataset, kinldy convert your dataset into this [Example](https://github.com/SCBB-LAB/Comparative-analysis-of-plant-TFBS-software/tree/main/Transcription_factor_motif_identification/CMF/cmfcode/ExampleDataSets) set.

#### 1.2 Usage
Run CMF program as shown below:
```
./cmf -w 7 -m 2 -d 1 -t 50 -i1 ExampleDataSets/youngOct4Bound.txt -i2 ExampleDataSets/youngOct4Control.txt -o seedsInfo.txt -f youngOct4Bound
```
**Final output:**
The output files are saved to the paraent diectory as  `output.txt`.

# Citation
If you use CMF in your research, please cite the following paper
[Co-regulation in embryonic stem cells via context-dependent binding of transcription factors](https://academic.oup.com/bioinformatics/article/29/17/2162/243962?login=true),


[Identification of Context-Dependent Motifs by Contrasting ChIP Binding Data](https://academic.oup.com/bioinformatics/article/26/22/2826/227908?login=true).
