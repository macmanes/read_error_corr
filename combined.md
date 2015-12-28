###Combined assembly of all 10 individuals

```
Trinity --seqType fq --max_memory 100G --trimmomatic \
--left SRR488153_1.fastq,SRR488152_1.fastq,SRR488154_1.fastq,SRR488155_1.fastq,SRR488332_1.fastq,SRR488333_1.fastq,SRR488334_1.fastq,SRR488335_1.fastq,SRR488337_1.fastq,SRR488338_1.fastq \
--right SRR488153_2.fastq,SRR488152_2.fastq,SRR488154_2.fastq,SRR488155_2.fastq,SRR488332_2.fastq,SRR488333_2.fastq,SRR488334_2.fastq,SRR488335_2.fastq,SRR488337_2.fastq,SRR488338_2.fastq \
--CPU 26 --output trinity_ALL_tuco --inchworm_cpu 10 --full_cleanup \
--quality_trimming_params "ILLUMINACLIP:/opt/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
```

### BUSCO evaluation

```
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -g trinity_ALL_tuco.Trinity.fasta -c 26 \
-m Trans -l /thomas2/matt_test/error/homo/vertebrata -o all_trinity_tuco
```

### Transrate evaluation

```
transrate -t 22 -a trinity_ALL_tuco.Trinity.fasta \
-l SRR488153_1.fastq,SRR488152_1.fastq,SRR488154_1.fastq,SRR488155_1.fastq,SRR488332_1.fastq,SRR488333_1.fastq,SRR488334_1.fastq,SRR488335_1.fastq,SRR488337_1.fastq,SRR488338_1.fastq \
-r SRR488153_2.fastq,SRR488152_2.fastq,SRR488154_2.fastq,SRR488155_2.fastq,SRR488332_2.fastq,SRR488333_2.fastq,SRR488334_2.fastq,SRR488335_2.fastq,SRR488337_2.fastq,SRR488338_2.fastq
```
