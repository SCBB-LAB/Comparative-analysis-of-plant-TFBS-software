# iFORM (Find Occurrence of Regulatory Motifs (iFORM)) 
This is a tool for scanning DNA sequences with TF motifs represented as PWMs, which integrates five traditional motif discovery programs, namely FIMO, Consensus, STORM, RSAT, and HOMER, through Fisher’s combined probability test.
To run iFORM on your input data, enter into the paraent directory:
```
cd iFORM
```

## 1.1 Data information
The software require the transcription factor specific Position Weight Matrix (PWM) and the genomic fasta file.
## 1.2 Usage
TO run iFORM successfully run following command:
```
./iForm --bgfile motif-file --o output_result Example/AT2G28810_pwm.txt Example/AT2G28810_pos.fa >> AT2G28810_out.txt
```
**Final output:** Final motif file is saved to `AT2G28810_out.txt` in the parent directory.
