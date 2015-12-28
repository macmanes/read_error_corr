###Assembly of a single individual

```bash
Trinity --seqType fq --max_memory 40G --trimmomatic \
--left SRR488338_1.fastq --right SRR488338_2.fastq \
--CPU 26 --output trinity_SRR488338_tuco --inchworm_cpu 6 --full_cleanup \
--quality_trimming_params "ILLUMINACLIP:/opt/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
```

###BUSCO eval

```bash
python3 /opt/BUSCO_v1.1b1/BUSCO_v1.1b1.py -g trinity_SRR488338_tuco.Trinity.fasta -c 26 \
-m Trans -l /thomas2/matt_test/error/homo/vertebrata -o trinity_tuco
```

###Transrate eval

```bash
transrate -t 26 -a trinity_SRR488338_tuco.Trinity.fasta \
-l SRR488338_1.fastq \
-r SRR488338_2.fastq
```
