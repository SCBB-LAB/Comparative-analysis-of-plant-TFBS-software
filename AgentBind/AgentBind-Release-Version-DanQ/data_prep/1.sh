fimo --oc example/ABF2_fimo example/ABF2_fimo/meme_m1.txt example/genomes/tair10/tair10.fa

#sed '$d' -i example/ABF2_fimo/fimo.tsv ; sed '$d' -i example/ABF2_fimo/fimo.tsv ; sed '$d' -i example/ABF2_fimo/fimo.tsv ; sed '$d' -i example/ABF2_fimo/fimo.tsv 

#python3 data_pre_software.py  --fimo_file example/ABF2_fimo/fimo.tsv --chipseq example/ABF2_narrow.bed --scope 1000 --resultdir example/ABF2 --datadir example --blockcore c --recorddir example/ABF2

#python2.7 ../model/train-transfer-learning.py  --data_dir example/ABF2/training --valid_dir example/ABF2/validation --n_train_samples 8974 --n_valid_samples 2238 --train_dir ABF2_out --seq_size 1000 --checkpoint_dir ABF2_out  --batch_size 128 --n_classes 2

