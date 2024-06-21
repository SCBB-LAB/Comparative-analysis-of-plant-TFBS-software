mkdir SRR8525028_chip
samtools view SRR8525028.bam | ./bin/MUSIC -preprocess SAM stdin SRR8525028_chip
samtools view SRR12022262.bam | ./bin/MUSIC -preprocess SAM stdin SRR12022262_input 
mkdir SRR8525028_chip/SRR8525028_sorted; mkdir SRR8525028_chip/SRR8525028_dedup

./bin/MUSIC -sort_reads SRR8525028_chip SRR8525028_chip/SRR8525028_sorted
./bin/MUSIC -sort_reads SRR12022262_input SRR12022262_input/sorted
./bin/MUSIC -remove_duplicates SRR8525028_chip/SRR8525028_sorted 2 SRR8525028_chip/SRR8525028_dedup
./bin/MUSIC -remove_duplicates SRR12022262_input/sorted 2 SRR12022262_input/dedup 

cd SRR8525028_chip

.././bin/MUSIC -get_multiscale_broad_ERs -chip ../SRR8525028_chip/SRR8525028_dedup -control ../SRR12022262_input -mapp Mappability_36bp -l_mapp 36 -begin_l 1000 -end_l 16000 -step 1.5 
cd ../
