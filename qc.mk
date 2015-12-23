#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#       ./qc.mk QC CPU=32 READ1=SRR797058_1.fastq READ2=SRR797058_2.fastq ASSEMBLY=1M.trinity_bfc31.Trinity.fasta RUN=1bfc31
#       ./qc.mk QC CPU=32 READ1=SRR797058_1.fastq READ2=SRR797058_2.fastq ASSEMBLY=5M.trinity_bfc31.Trinity.fasta RUN=5bfc31
#       ./qc.mk QC CPU=32 READ1=SRR797058_1.fastq READ2=SRR797058_2.fastq ASSEMBLY=10M.trinity_bfc31.Trinity.fasta RUN=10bfc31
#       ./qc.mk QC CPU=32 READ1=SRR797058_1.fastq READ2=SRR797058_2.fastq ASSEMBLY=20M.trinity_bfc31.Trinity.fasta RUN=20bfc31
#       ./qc.mk QC CPU=32 READ1=SRR797058_1.fastq READ2=SRR797058_2.fastq ASSEMBLY=60M.trinity_rcorr31.Trinity.fasta RUN=60rcorr
#       ./qc.mk QC CPU=32 READ1=SRR797058_1.fastq READ2=SRR797058_2.fastq ASSEMBLY=80M.trinity_rcorr31.Trinity.fasta RUN=80rcorr


MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
READ1=SRR797058_1.fastq
READ2=SRR797058_2.fastq
ASSEMBLY=Trinity.fasta
RUN=
BUSCO ?= ${shell which BUSCO_v1.1b1.py}
BUSCODIR := $(dir $(firstword $(BUSCO)))
SAMP=20

all:QC diginorm
.PHONY:QC diginorm
.DELETE_ON_ERROR:


download_reads:
	mkdir -p ${DIR}/reads
	cd ${DIR}/reads && \
	curl -L https://s3.amazonaws.com/macmanes_general/SRR797058_2.fastq.gz | gzip -cd > SRR797058_2.fastq && \
	curl -L https://s3.amazonaws.com/macmanes_general/SRR797058_1.fastq.gz | gzip -cd > SRR797058_1.fastq

QC:
	mkdir -p ${DIR}/QC
	cd ${DIR}/QC && \
	python3 ${BUSCODIR}BUSCO_v1.1b1.py -o RAW_BUSCO_${ASSEMBLY} -g ../assemblies/${ASSEMBLY} -m Trans --cpu $(CPU) -l ${BUSCODIR}/vertebrata && \
	transrate --output RAW_TRANSRATE_${ASSEMBLY} -a ../assemblies/${ASSEMBLY} --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	python3 ${BUSCODIR}BUSCO_v1.1b1.py -o GOOD_BUSCO_${ASSEMBLY} -g RAW_BUSCO_${ASSEMBLY}/good*fasta -m Trans --cpu $(CPU) -l ${BUSCODIR}/vertebrata && \
	kallisto index -i kallisto.idx ../assemblies/${ASSEMBLY} && \
	kallisto quant -t $(CPU) -i kallisto.idx -o kallisto_orig -b 100 ${DIR}/reads/${READ1} ${DIR}/reads/${READ2} && \
	~/salmon-0.5.1/bin/salmon index -t ../assemblies/${ASSEMBLY} -i salmon.idx --type quasi -k 31 && \
	~/salmon-0.5.1/bin/salmon quant -p $(CPU) -i transcripts2_index -l MSR -1 ${DIR}/reads/${READ1} -2 ${DIR}/reads/${READ2} -o salmon_orig && \
	awk '1>$5{next}1' kallisto_orig/abundance.tsv | awk '{print $1}' > list && \
	awk '1>$3{next}1' salmon_orig/quant.sf | sed  '1,10d' | awk '{print $1}' > list2 && \
	sed -i ':begin;$!N;/[ACTGNn-]\n[ACTGNn-]/s/\n//;tbegin;P;D' ../assemblies/${ASSEMBLY} && \
	for i in $(cat list_final); do grep --no-group-separator --max-count=1 -A1 -w $i ../assemblies/${ASSEMBLY} >> ${RUN}.Trinity.fasta; done && \
	python3 ${BUSCODIR}BUSCO_v1.1b1.py -o TPM_FILT_BUSCO_${RUN} -g ${RUN}.Trinity.fasta -m Trans --cpu $(CPU) -l ${BUSCODIR}/vertebrata && \
	transrate -o ${SAMP}M.bless55 --output TPM_FILT_TRANSRATE_${RUN} -a ${RUN}.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU)

diginorm:
	mkdir -p ${DIR}/diginorm
	cd -p ${DIR}/diginorm && \
	seqtk sample -s4102340 ${READ1} ${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}a.subsamp_1.fastq && \
	seqtk sample -s4102340 ${READ2} ${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}a.subsamp_2.fastq && \
	seqtk mergepe 20a.subsamp_1.fastq 20a.subsamp_2.fastq > 20a.interleave.fastq && \
	normalize-by-median.py --max-memory-usage 4e9 -C 30 -o 20a.norm.fastq 20a.interleave.fastq && \
	split-paired-reads.py 20a.norm.fastq && \
	mv 20a.norm.fastq.1 20a.norm.1.fastq && \
	mv 20a.norm.fastq.2 20a.norm.2.fastq && \
	seqtk sample -s222323440 ${READ1} ${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}b.subsamp_1.fastq && \
	seqtk sample -s222323440 ${READ2} ${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}b.subsamp_2.fastq && \
	seqtk mergepe 20b.subsamp_1.fastq 20b.subsamp_2.fastq > 20b.interleave.fastq && \
	normalize-by-median.py --max-memory-usage 4e9 -C 30 -o 20b.norm.fastq 20b.interleave.fastq && \
	split-paired-reads.py 20b.norm.fastq && \
	mv 20b.norm.fastq.1 20b.norm.1.fastq && \
	mv 20b.norm.fastq.2 20b.norm.2.fastq && \
	seqtk sample -s76232344098 ${READ1} ${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}c.subsamp_1.fastq && \
	seqtk sample -s76232344098 ${READ2} ${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}c.subsamp_2.fastq && \
	seqtk mergepe 20c.subsamp_1.fastq 20c.subsamp_2.fastq > 20c.interleave.fastq && \
	normalize-by-median.py --max-memory-usage 4e9 -C 30 -o 20c.norm.fastq 20c.interleave.fastq && \
	split-paired-reads.py 20c.norm.fastq && \
	mv 20c.norm.fastq.1 20c.norm.1.fastq && \
	mv 20c.norm.fastq.2 20c.norm.2.fastq && \
	Trinity --seqType fq --max_memory 20G --trimmomatic --left 20a.norm.1.fastq --right 20a.norm.2.fastq --CPU $(CPU) --output ${SAMP}mus_A_C30_trinity --inchworm_cpu 6 --full_cleanup --quality_trimming_params "ILLUMINACLIP:/home/ubuntu/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --seqType fq --max_memory 20G --trimmomatic --left 20b.norm.1.fastq --right 20b.norm.2.fastq --CPU $(CPU) --output ${SAMP}mus_B_C30_trinity --inchworm_cpu 6 --full_cleanup --quality_trimming_params "ILLUMINACLIP:/home/ubuntu/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --seqType fq --max_memory 20G --trimmomatic --left 20c.norm.1.fastq --right 20c.norm.2.fastq --CPU $(CPU) --output ${SAMP}mus_C_C30_trinity --inchworm_cpu 6 --full_cleanup --quality_trimming_params "ILLUMINACLIP:/home/ubuntu/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	python3 ${BUSCODIR}BUSCO_v1.1b1.py  -g ${SAMP}mus_A_C30_trinity.Trinity.fasta -m Trans --cpu $(CPU) -l ${BUSCODIR}/vertebrata && \
	python3 ${BUSCODIR}BUSCO_v1.1b1.py  -g ${SAMP}mus_B_C30_trinity.Trinity.fasta -m Trans --cpu $(CPU) -l ${BUSCODIR}/vertebrata && \
	python3 ${BUSCODIR}BUSCO_v1.1b1.py  -g ${SAMP}mus_C_C30_trinity.Trinity.fasta -m Trans --cpu $(CPU) -l ${BUSCODIR}/vertebrata && \
	transrate --output ${SAMP}A_TRANSRATE_${ASSEMBLY} -a ${SAMP}mus_A_C30_trinity.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	transrate --output ${SAMP}B_TRANSRATE_${ASSEMBLY} -a ${SAMP}mus_B_C30_trinity.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	transrate --output ${SAMP}C_TRANSRATE_${ASSEMBLY} -a ${SAMP}mus_C_C30_trinity.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU)
