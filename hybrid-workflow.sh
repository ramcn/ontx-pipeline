cd dataset
flappie  ../fast5_files/ > ont.fasta
spades.py --careful --pe1-1 Kp2146_paired_1.fastq.gz --pe1-2 Kp2146_paired_2.fastq.gz -o spades -t 16
darwin-xl ecoli_spades_50.fasta reads_50.fasta 0 0  > out.sam
nanopolish variants --consensus -o polished.vcf -w "tig00000001:200000-202000" -r reads.fasta -b reads.sorted.bam -g draft.fa 
