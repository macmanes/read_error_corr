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
	python3 ${BUSCODIR}BUSCO_v1.1b1.py -o GOOD_BUSCO_${ASSEMBLY} -g RAW_TRANSRATE_${ASSEMBLY}/*/good*fasta -m Trans --cpu $(CPU) -l ${BUSCODIR}/vertebrata && \
	kallisto index -i kallisto.idx ../assemblies/${ASSEMBLY} && \
	kallisto quant -t $(CPU) -i kallisto.idx -o kallisto_orig -b 100 ${DIR}/reads/${READ1} ${DIR}/reads/${READ2} && \
	salmon index -t ../assemblies/${ASSEMBLY} -i salmon.idx --type quasi -k 31 && \
	salmon quant -p $(CPU) -i salmon.idx -l MSR -1 ${DIR}/reads/${READ1} -2 ${DIR}/reads/${READ2} -o salmon_orig && \
	awk '1>$$5{next}1' kallisto_orig/abundance.tsv | awk '{print $$1}' > list && \
	awk '1>$$3{next}1' salmon_orig/quant.sf | sed  '1,10d' | awk '{print $$1}' > list2 && \
	sed -i ':begin;$!N;/[ACTGNn-]\n[ACTGNn-]/s/\n//;tbegin;P;D' ../assemblies/${ASSEMBLY} && \
	cat list list2 | sort -u > list_final && \
	for i in $$(cat list_final); do grep --no-group-separator --max-count=1 -A1 -w $$i ../assemblies/${ASSEMBLY} >> ${RUN}.Trinity.fasta; done && \
	python3 ${BUSCODIR}BUSCO_v1.1b1.py -o TPM_FILT_BUSCO_${RUN} -g ${RUN}.Trinity.fasta -m Trans --cpu $(CPU) -l ${BUSCODIR}/vertebrata && \
	transrate --output TPM_FILT_TRANSRATE_${RUN} -a ${RUN}.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU)

diginorm:
	mkdir -p ${DIR}/diginorm
	cd ${DIR}/diginorm && \
	seqtk sample -s4102340 ../reads/${READ1} ${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}a.subsamp_1.fastq && \
	seqtk sample -s4102340 ../reads/${READ2} ${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}a.subsamp_2.fastq && \
	seqtk mergepe ${SAMP}a.subsamp_1.fastq ${SAMP}a.subsamp_2.fastq > ${SAMP}a.interleave.fastq && \
	normalize-by-median.py --max-memory-usage 4e9 -C 30 -o ${SAMP}a.norm.fastq ${SAMP}a.interleave.fastq && \
	split-paired-reads.py ${SAMP}a.norm.fastq && \
	mv ${SAMP}a.norm.fastq.1 ${SAMP}a.norm.1.fastq && \
	mv ${SAMP}a.norm.fastq.2 ${SAMP}a.norm.2.fastq && \
	seqtk sample -s222323440 ${READ1} ../reads/${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}b.subsamp_1.fastq && \
	seqtk sample -s222323440 ${READ2} ../reads/${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}b.subsamp_2.fastq && \
	seqtk mergepe ${SAMP}b.subsamp_1.fastq ${SAMP}b.subsamp_2.fastq > ${SAMP}b.interleave.fastq && \
	normalize-by-median.py --max-memory-usage 4e9 -C 30 -o ${SAMP}b.norm.fastq ${SAMP}b.interleave.fastq && \
	split-paired-reads.py ${SAMP}b.norm.fastq && \
	mv ${SAMP}b.norm.fastq.1 ${SAMP}b.norm.1.fastq && \
	mv ${SAMP}b.norm.fastq.2 ${SAMP}b.norm.2.fastq && \
	seqtk sample -s76232344098 ${READ1} ../reads/${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}c.subsamp_1.fastq && \
	seqtk sample -s76232344098 ${READ2} ../reads/${SAMP}000000 | sed 's_ H_-H_g' > ${SAMP}c.subsamp_2.fastq && \
	seqtk mergepe ${SAMP}c.subsamp_1.fastq ${SAMP}c.subsamp_2.fastq > ${SAMP}c.interleave.fastq && \
	normalize-by-median.py --max-memory-usage 4e9 -C 30 -o ${SAMP}c.norm.fastq ${SAMP}c.interleave.fastq && \
	split-paired-reads.py ${SAMP}c.norm.fastq && \
	mv ${SAMP}c.norm.fastq.1 ${SAMP}c.norm.1.fastq && \
	mv ${SAMP}c.norm.fastq.2 ${SAMP}c.norm.2.fastq && \
	Trinity --seqType fq --max_memory 20G --trimmomatic --left ${SAMP}a.norm.1.fastq --right ${SAMP}a.norm.2.fastq --CPU $(CPU) --output ${SAMP}mus_A_C30_trinity --inchworm_cpu 6 --full_cleanup --quality_trimming_params "ILLUMINACLIP:/home/ubuntu/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --seqType fq --max_memory 20G --trimmomatic --left ${SAMP}b.norm.1.fastq --right ${SAMP}b.norm.2.fastq --CPU $(CPU) --output ${SAMP}mus_B_C30_trinity --inchworm_cpu 6 --full_cleanup --quality_trimming_params "ILLUMINACLIP:/home/ubuntu/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --seqType fq --max_memory 20G --trimmomatic --left ${SAMP}c.norm.1.fastq --right ${SAMP}c.norm.2.fastq --CPU $(CPU) --output ${SAMP}mus_C_C30_trinity --inchworm_cpu 6 --full_cleanup --quality_trimming_params "ILLUMINACLIP:/home/ubuntu/trinityrnaseq/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	python3 ${BUSCODIR}BUSCO_v1.1b1.py  -g ${SAMP}mus_A_C30_trinity.Trinity.fasta -m Trans --cpu $(CPU) -l ${BUSCODIR}/vertebrata && \
	python3 ${BUSCODIR}BUSCO_v1.1b1.py  -g ${SAMP}mus_B_C30_trinity.Trinity.fasta -m Trans --cpu $(CPU) -l ${BUSCODIR}/vertebrata && \
	python3 ${BUSCODIR}BUSCO_v1.1b1.py  -g ${SAMP}mus_C_C30_trinity.Trinity.fasta -m Trans --cpu $(CPU) -l ${BUSCODIR}/vertebrata && \
	transrate --output ${SAMP}A_TRANSRATE_${ASSEMBLY} -a ${SAMP}mus_A_C30_trinity.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	transrate --output ${SAMP}B_TRANSRATE_${ASSEMBLY} -a ${SAMP}mus_B_C30_trinity.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	transrate --output ${SAMP}C_TRANSRATE_${ASSEMBLY} -a ${SAMP}mus_C_C30_trinity.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU)
