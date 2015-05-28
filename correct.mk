#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
SAMP=10


all: scripts download_reads subsamp_reads reference raw lighter bless sga bfc seecer


scripts:
	@echo Downloading Scripts
	mkdir -p ${DIR}/scripts
	cd ${DIR}/scripts && \
	wget https://raw.githubusercontent.com/macmanes/trimming_paper/master/scripts/subsampler.py


download_reads:
	mkdir -p ${DIR}/reads
	cd ${DIR}/reads && \
	gzip -cd <(wget -qO- ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR797/SRR797058/SRR797058_1.fastq.gz) > SRR797058_1.fastq && \
	gzip -cd <(wget -qO- ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR797/SRR797058/SRR797058_2.fastq.gz) > SRR797058_2.fastq
subsamp_reads:
	cd ${DIR}/reads && \
	python ${DIR}/scripts/subsampler.py $(SAMP)000000 SRR797058_1.fastq SRR797058_2.fastq && \
	sed -i 's_ H_-H_g' subsamp_{1,2}.fastq

reference:
	mkdir -p ${DIR}/genome
	cd ${DIR}/genome && \
	wget ftp://ftp.ensembl.org/pub/release-79/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.1.fa.gz && \
	gzip -d Mus_musculus.GRCm38.dna.chromosome.1.fa.gz && \
	bwa index -p mus Mus_musculus.GRCm38.dna.chromosome.1.fa

raw:
	mkdir -p ${DIR}/raw
	cd ${DIR}/raw && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq > $(SAMP)M.raw.sam

lighter:
	mkdir -p ${DIR}/lighter$(SAMP)M
	cd ${DIR}/lighter$(SAMP)M && \
	lighter -K 31 60000000 -r ${DIR}/reads/subsamp_1.fastq -r ${DIR}/reads/subsamp_2.fastq -t $(CPU) && \
	bwa mem -t $(CPU) ${DIR}/genome/mus subsamp_1.cor.fq subsamp_2.cor.fq > $(SAMP)M.lighter.sam



bless:
	mkdir -p ${DIR}/bless$(SAMP)M
	mkdir -p ${DIR}/bless$(SAMP)M/kmc/bin/
	cp /home/ubuntu/v0p24/kmc/bin/kmc /mnt/bless$(SAMP)M/kmc/bin/kmc
	cd ${DIR}/bless$(SAMP)M && \
	mpirun -np $(CPU) bless -read1 ${DIR}/reads/subsamp_1.fastq -read2 ${DIR}/reads/subsamp_2.fastq -prefix $(SAMP)M_bless55 -kmerlength 55 && \
	mpirun -np $(CPU) bless -read1 ${DIR}/reads/subsamp_1.fastq -read2 ${DIR}/reads/subsamp_2.fastq -prefix $(SAMP)M_bless33 -kmerlength 33 && \
	bwa mem -t $(CPU) ${DIR}/genome/mus $(SAMP)M_bless55.1.corrected.fastq $(SAMP)M_bless55.2.corrected.fastq > $(SAMP)M.bless55.sam && \
	bwa mem -t $(CPU) ${DIR}/genome/mus $(SAMP)M_bless33.1.corrected.fastq $(SAMP)M_bless33.2.corrected.fastq > $(SAMP)M.bless33.sam


sga:
	mkdir -p ${DIR}/sga$(SAMP)M
	cd ${DIR}/sga$(SAMP)M && \
	sga preprocess -p 1 ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq | gzip -1 > out.pe.fq.gz && \
	sga index -a ropebwt -t 8 --no-reverse out.pe.fq.gz && \
	sga correct -t $(CPU) -k 55 --learn out.pe.fq.gz -o sga.55.fq && \
	sga correct -t $(CPU) -k 33 --learn out.pe.fq.gz -o sga.33.fq && \
	split-paired-reads.py sga.55.fq && \
	split-paired-reads.py sga.33.fq && \
	bwa mem -t $(CPU) ${DIR}/genome/mus sga.55.fq.1 sga.55.fq.2 > $(SAMP)M.sga55.sam && \
	bwa mem -t $(CPU) ${DIR}/genome/mus sga.33.fq.1 sga.33.fq.2 > $(SAMP)M.sga33.sam


bfc:
	mkdir -p ${DIR}/bfc$(SAMP)M
	cd ${DIR}/bfc$(SAMP)M && \
	interleave-reads.py ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq -o inter.fq && \
	bfc -s 50m -k55 -t 16 inter.fq > bfc55.corr.fq && \
	bfc -s 50m -k33 -t 16 inter.fq > bfc33.corr.fq && \
	split-paired-reads.py bfc55.corr.fq && \
	split-paired-reads.py bfc33.corr.fq && \
	bwa mem -t $(CPU) ${DIR}/genome/mus bfc55.corr.fq.1 bfc55.corr.fq.2 > $(SAMP)M.bfc55.sam && \
	bwa mem -t $(CPU) ${DIR}/genome/mus bfc33.corr.fq.1 bfc33.corr.fq.2 > $(SAMP)M.bfc33.sam

seecer:
	mkdir -p ${DIR}/seecer$(SAMP)M
	cd ${DIR}/seecer$(SAMP)M && \
	run_seecer.sh -t . -k 31 ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq
