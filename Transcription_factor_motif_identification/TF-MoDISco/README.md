# TF-MoDISco
TF-MoDISco (Transcription Factor MotifDiscovery from Importance Scores), a novel algorithm that leverages per-base importance scores to simultaneously incorporate information from all neurons in the network and generate high-quality, consolidated, non-redundant motifs.

## Python programs for predicting Transcription factor sites (TFBSs) or motif identification and performing DeepLIFT.
## Dependencies
The program requires:
  * python==3.7.13
  * tensorflow-gpu==2.0.0 (for model training)
  * tensorfow==1.14.0 (for DeepLIFT, TF-MoDISco and model predicting)
  * deeplift
  * modisco
  * keras==2.3.1.0.2
  * scikit-learn==1.0.2
  * pandas 
  * numpy 
  * the [bedtools](https://bedtools.readthedocs.io/en/latest/) software

## Tutorial
###  Usage
```
python3 modisco_test.py <input fasta file> <species> <tf>
```

It should be noted ```<species>``` that one is chosen from 'Zea_mays_models','Arabidopsis_models' and 'Oryza_sativa_models'.
It should be noted``` <tf>``` that one is chosen from the tf names of selected species.
After running the program, a dir about tf-modisco results will be generated in the current folder.
We here provide a `ATHB25_pos.fa` file and employed one of models of Arabidopsis thaliana for an example: 

```
python3 modisco2.py Example/ATHB25_pos.fa models/Arabidopsis_models ATHB25
```
```
### Citation
* Huang, G., et al. Densely Connected Convolutional Networks. IEEE Computer Society 2016.
* Shrikumar, A., Greenside, P. and Kundaje, A. Learning Important Features Through Propagating Activation Differences. 2017.
* Shrikumar, A., et al. Technical Note on Transcription Factor Motif Discovery from Importance Scores (TF-MoDISco) version 0.5.6.5. In.; 2018. p. arXiv:1811.00416.
