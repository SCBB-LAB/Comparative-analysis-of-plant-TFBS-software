#python3.8 deepRAM.py --test_data Example/$1_B.seq --data_type DNA --predict_only True --model_path Example/$1.pkl --motif True --motif_dir motifs --out_file Example/$1_prediction.txt  --Embedding False --Conv True --RNN False --conv_layers 1

ls motifs/pwm*.txt | sed 's+.txt++' | while read i ; do python3.8 x6.py $i Example/$1_pos.fa Example/$1_full ; done

cp motifs/* Example/pwm_logo/$1/.
