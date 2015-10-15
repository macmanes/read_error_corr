#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	for i in 10 20 30 40 50 60 70 80 90 100; do ./correct.mk main SAMP=$i CPU=36; done
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
SAMP=10
READ1=SRR797058_1.fastq.gz
READ2=SRR797058_2.fastq.gz


prep: download_reads scripts reference
main: subsamp_reads raw trinity_raw lighter lighter_trinity bless bless_trinity sga sga_trinity bfc bfc_trinity seecer seecer_trinity \
	rcorrector 

.DELETE_ON_ERROR:

scripts:
	@echo Downloading Scripts
	mkdir -p ${DIR}/scripts
	cd ${DIR}/scripts && \
	wget https://raw.githubusercontent.com/macmanes/read_error_corr/master/barcodes.fa

download_reads:
	mkdir -p ${DIR}/reads
	cd ${DIR}/reads && \
	curl -LO ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR797/SRR797058/SRR797058_1.fastq.gz && \
	curl -LO ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR797/SRR797058/SRR797058_2.fastq.gz

reference:
	mkdir -p ${DIR}/genome
	cd ${DIR}/genome && \
	wget ftp://ftp.ensembl.org/pub/release-79/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.1.fa.gz && \
	bwa index -p mus Mus_musculus.GRCm38.dna.chromosome.1.fa.gz

subsamp_reads:
	cd ${DIR}/reads && \
	seqtk sample -s102340 ${READ1} ${SAMP}000000 > subsamp_1.fastq && \
	seqtk sample -s102340 ${READ2} ${SAMP}000000 > subsamp_2.fastq && \
	sed -i 's_ H_-H_g' subsamp_{1,2}.fastq

raw:${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq
	mkdir -p ${DIR}/raw
	cd ${DIR}/raw && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq > ${SAMP}M.raw.sam

lighter:${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq
	mkdir -p ${DIR}/lighter${SAMP}M
	cd ${DIR}/lighter${SAMP}M && \
	lighter -K 31 60000000 -r ${DIR}/reads/subsamp_1.fastq -r ${DIR}/reads/subsamp_2.fastq -t $(CPU) && \
	bwa mem -t $(CPU) ${DIR}/genome/mus subsamp_1.cor.fq subsamp_2.cor.fq > ${SAMP}M.lighter.sam && \
	k8 ~/bfc/errstat.js ${DIR}/lighter${SAMP}M/${SAMP}M.lighter.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.lighter.out

lighter_trinity:
	cd ${DIR}/lighter${SAMP}M && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left subsamp_1.cor.fq --right subsamp_2.cor.fq --CPU $(CPU) --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"

bless:
	mkdir -p ${DIR}/bless${SAMP}M
	cd ${DIR}/bless${SAMP}M && \
	mpirun -np $(CPU) bless -notrim -read1 ${DIR}/reads/subsamp_1.fastq -read2 ${DIR}/reads/subsamp_2.fastq -prefix ${SAMP}M_bless55 -kmerlength 55 && \
	mpirun -np $(CPU) bless -notrim -read1 ${DIR}/reads/subsamp_1.fastq -read2 ${DIR}/reads/subsamp_2.fastq -prefix ${SAMP}M_bless31 -kmerlength 31 && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${SAMP}M_bless55.1.corrected.fastq ${SAMP}M_bless55.2.corrected.fastq > ${SAMP}M.bless55.sam && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${SAMP}M_bless31.1.corrected.fastq ${SAMP}M_bless31.2.corrected.fastq > ${SAMP}M.bless31.sam && \
	k8 ~/bfc/errstat.js ${DIR}/bless${SAMP}M/${SAMP}M.bless55.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bless55.out && \
	k8 ~/bfc/errstat.js ${DIR}/bless${SAMP}M/${SAMP}M.bless31.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bless31.out

bless_trinity:
	cd ${DIR}/bless${SAMP}M && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left ${SAMP}M_bless55.1.corrected.fastq --right ${SAMP}M_bless55.2.corrected.fastq --CPU $(CPU) --output trinity_bless55 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left ${SAMP}M_bless31.1.corrected.fastq --right ${SAMP}M_bless31.2.corrected.fastq --CPU $(CPU) --output trinity_bless31 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"

sga:
	mkdir -p ${DIR}/sga${SAMP}M
	cd ${DIR}/sga${SAMP}M && \
	sga preprocess -p 1 ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq | gzip -1 > out.pe.fq.gz && \
	sga index -a ropebwt -t 8 --no-reverse out.pe.fq.gz && \
	sga correct -t $(CPU) -k 55 --learn out.pe.fq.gz -o sga.55.fq && \
	sga correct -t $(CPU) -k 31 --learn out.pe.fq.gz -o sga.31.fq && \
	bwa mem -p -t $(CPU) ${DIR}/genome/mus sga.55.fq > ${SAMP}M.sga55.sam && \
	bwa mem -p -t $(CPU) ${DIR}/genome/mus sga.31.fq > ${SAMP}M.sga31.sam && \
	k8 ~/bfc/errstat.js ${DIR}/sga${SAMP}M/${SAMP}M.sga55.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.sga55.out && \
	k8 ~/bfc/errstat.js ${DIR}/sga${SAMP}M/${SAMP}M.sga31.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.sga31.out

sga_trinity:
	cd ${DIR}/sga${SAMP}M && \
	split-paired-reads.py ${DIR}/genome/mus sga.55.fq && \
	split-paired-reads.py ${DIR}/genome/mus sga.31.fq && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left sga.55.fq.1 --right sga.55.fq.2 --CPU $(CPU) --output trinity_sga55 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left sga.31.fq.1 --right sga.31.fq.2 --CPU $(CPU) --output trinity_sga31 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"


bfc:
	mkdir -p ${DIR}/bfc${SAMP}M
	cd ${DIR}/bfc${SAMP}M && \
	seqtk mergepe ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq > inter.fq && \
	bfc -s 50m -k55 -t $(CPU) inter.fq > bfc55.corr.fq && \
	bfc -s 50m -k31 -t $(CPU) inter.fq > bfc31.corr.fq && \
	bwa mem -p -t $(CPU) ${DIR}/genome/mus bfc55.corr.fq > ${SAMP}M.bfc55.sam && \
	bwa mem -p -t $(CPU) ${DIR}/genome/mus bfc31.corr.fq > ${SAMP}M.bfc31.sam && \
	k8 ~/bfc/errstat.js ${DIR}/bfc${SAMP}M/${SAMP}M.bfc55.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bfc55.out && \
	k8 ~/bfc/errstat.js ${DIR}/bfc${SAMP}M/${SAMP}M.bfc31.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bfc31.out

bfc_trinity:
	cd ${DIR}/bfc${SAMP}M && \
	split-paired-reads.py ${DIR}/genome/mus bfc55.corr.fq && \
	split-paired-reads.py ${DIR}/genome/mus bfc31.corr.fq && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left bfc55.corr.fq.1 --right bfc55.corr.fq.2 --CPU $(CPU) --output trinity_bfc55 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left bfc31.corr.fq.1 --right bfc31.corr.fq.2 --CPU $(CPU) --output trinity_bfc31 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"


seecer:
	mkdir -p ${DIR}/seecer${SAMP}M
	cd ${DIR}/seecer${SAMP}M && \
	run_seecer.sh -t . -k 31 ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${DIR}/reads/subsamp_1.fastq_corrected.fa ${DIR}/reads/subsamp_2.fastq_corrected.fa > ${SAMP}M.seecer.sam && \
	k8 ~/bfc/errstat.js ${DIR}/seecer${SAMP}M/${SAMP}M.seecer.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.seecer.out

seecer_trinity:
	cd ${DIR}/seecer${SAMP}M && \
	Trinity --seqType fa --max_memory 50G --trimmomatic --left ${DIR}/reads/subsamp_1.fastq_corrected.fa --right ${DIR}/reads/subsamp_2.fastq_corrected.fa --CPU $(CPU) --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"

rcorrector:
	mkdir -p ${DIR}/rcorr${SAMP}M
	cd ${DIR}/rcorr${SAMP}M && \
	perl /share/Rcorrector/run_rcorrector.pl -t $(CPU) -k 25 -1 ${DIR}/reads/subsamp_1.fastq -2 ${DIR}/reads/subsamp_2.fastq


trinity_raw:${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq
	mkdir -p ${DIR}/trinity_${SAMP}M
	cd ${DIR}/trinity_${SAMP}M && \
	Trinity --seqType fq --max_memory 50G --trimmomatic \
	--left $< \
	--right $(word 2,$^) \
	--CPU $(CPU) --output trinity_${SAMP}M.P2.raw --inchworm_cpu 10 --full_cleanup \
	--quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
