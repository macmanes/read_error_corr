###Donload and eval original Oreochromis niloticus assembly

```bash

curl -LO http://www.ncbi.nlm.nih.gov/Traces/wgs/?download=GAID01.1.fsa_nt.gz
mv ?download=GAID01.1.fsa_nt.gz GAID01.1.fsa_nt.gz
gzip -d GAID01.1.fsa_nt.gz

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR797/SRR797490/SRR797490_1.fastq.gz \
| gzip -d | sed 's_ __g' > SRR797490_1.fastq

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR797/SRR797490/SRR797490_2.fastq.gz \
| gzip -d | sed 's_ __g' > SRR797490_2.fastq

python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -g GAID01.1.fsa_nt -c 32 -o RAW_GAID01 -m Trans -l ~/BUSCO_v1.1b1/vertebrata

transrate -t 32 -a GAID01.1.fsa_nt -o raw_GAID01 \
-l SRR797490_1.fastq \
-r SRR797490_2.fastq
```


###RCorrector reads

```bash
perl run_rcorrector.pl -k 31 -t 32 \
-1 SRR797490_1.fastq \
-2 SRR797490_2.fastq
```

###Trinity Oreochromis niloticus

```bash
Trinity --seqType fq --max_memory 20G --trimmomatic --CPU 30 \
--left SRR797490_1.cor.fq --right SRR797490_2.cor.fq --output trinity_GAID01_P2_Rcorr \
--quality_trimming_params "ILLUMINACLIP:/home/ubuntu/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
```

###Eval Trinity transcriptome Oreochromis niloticus


```bash
python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -g trinity_GAID01_P2_Rcorr.Trinity.fasta \
-c 32 -o trinity -m Trans -l ~/BUSCO_v1.1b1/vertebrata


transrate -t 32 -a trinity_GAID01_P2_Rcorr.Trinity.fasta -o fish_trinity \
-l SRR797490_1.cor.fq \
-r SRR797490_2.cor.fq
```

###Eval good assembly from Trinity

```bash
python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -g fish_trinity/*/good*fasta \
-c 32 -o fish_good -m Trans -l ~/BUSCO_v1.1b1/vertebrata
```

### TPM Filter Trinity Oreochromis assembly

```bash
kallisto index -i kallisto.idx trinity_GAID01_P2_Rcorr.Trinity.fasta
kallisto quant -t 32 -i kallisto.idx -o kallisto_orig -b 100 SRR797490_1.cor.fq SRR797490_2.cor.fq

~/salmon-0.5.1/bin/salmon index -t trinity_GAID01_P2_Rcorr.Trinity.fasta -i salmon.idx --type quasi -k 31
~/salmon-0.5.1/bin/salmon quant -p 32 -i salmon.idx -l MSR -1 SRR797490_1.cor.fq -2 SRR797490_2.cor.fq -o salmon_orig

awk '1>$5{next}1' kallisto_orig/abundance.tsv | awk '{print $1}' > kallist
awk '1>$3{next}1' salmon_orig/quant.sf | sed  '1,10d' | awk '{print $1}' > salist
cat kallist salist | sort -u > uniq_list
sed -i ':begin;N;/[ACTGNn-]\n[ACTGNn-]/s/\n//;tbegin;P;D' trinity_GAID01_P2_Rcorr.Trinity.fasta

for i in $(cat uniq_list);
   do grep --no-group-separator --max-count=1 -A1 -w $i trinity_GAID01_P2_Rcorr.Trinity.fasta >> highexp_GAID01_P2_Rcorr.Trinity.fasta;
done

```

###Eval TPM filter Trinity Oreochromis assembly

```bash
python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -g highexp_GAID01_P2_Rcorr.Trinity.fasta \
-c 32 -o highexp -m Trans -l ~/BUSCO_v1.1b1/vertebrata


transrate -t 32 -a highexp_GAID01_P2_Rcorr.Trinity.fasta -o fish_highexp \
-l SRR797490_1.cor.fq \
-r SRR797490_2.cor.fq
```
### TPM Filter Transrate GOOD assembly

```bash
kallisto index -i kallisto.idx fish_trinity/trinity_GAID01_P2_Rcorr.Trinity/good.trinity_GAID01_P2_Rcorr.Trinity.fasta
kallisto quant -t 32 -i kallisto.idx -o kallisto_orig -b 100 SRR797490_1.cor.fq SRR797490_2.cor.fq

~/salmon-0.5.1/bin/salmon index -t fish_trinity/trinity_GAID01_P2_Rcorr.Trinity/good.trinity_GAID01_P2_Rcorr.Trinity.fasta -i salmon.idx --type quasi -k 31
~/salmon-0.5.1/bin/salmon quant -p 32 -i salmon.idx -l MSR -1 SRR797490_1.cor.fq -2 SRR797490_2.cor.fq -o salmon_orig

awk '1>$5{next}1' kallisto_orig/abundance.tsv | awk '{print $1}' > kallist
awk '1>$3{next}1' salmon_orig/quant.sf | sed  '1,10d' | awk '{print $1}' > salist
cat kallist salist | sort -u > uniq_list
sed -i ':begin;N;/[ACTGNn-]\n[ACTGNn-]/s/\n//;tbegin;P;D' fish_trinity/trinity_GAID01_P2_Rcorr.Trinity/good.trinity_GAID01_P2_Rcorr.Trinity.fasta

for i in $(cat uniq_list);
   do grep --no-group-separator --max-count=1 -A1 -w $i fish_trinity/trinity_GAID01_P2_Rcorr.Trinity/good.trinity_GAID01_P2_Rcorr.Trinity.fasta >> highexp_GOOD_GAID01_P2_Rcorr.Trinity.fasta;
done

```

###Eval TPM Filter Transrate GOOD assembly

```bash
python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -g highexp_GOOD_GAID01_P2_Rcorr.Trinity.fasta \
-c 32 -o highexp_good -m Trans -l ~/BUSCO_v1.1b1/vertebrata


transrate -t 32 -a highexp_GOOD_GAID01_P2_Rcorr.Trinity.fasta -o fish_highexp_good \
-l SRR797490_1.cor.fq \
-r SRR797490_2.cor.fq
```
