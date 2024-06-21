ls *_sorted.bam | while read i ; do samtools index $i ; done
ls *_sorted.bam | while read i ; do samtools flagstat $i ; done
ls *_sorted.bam | sed 's+_sorted.bam++' | while read i ; do bamCoverage --bam $i\_sorted.bam -o $i ; done
