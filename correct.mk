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
BFC ?= ${shell which bfc}
BFCDIR := $(dir $(firstword $(BFC)))
BLESS ?= ${shell which bless}
KMC := $(dir $(firstword $(BFC)))
READ1=SRR797058_1.fastq.gz
READ2=SRR797058_2.fastq.gz


prep: download_reads scripts reference
main: setup subsamp_reads raw trinity_raw lighter lighter_trinity bless bless_trinity sga sga_trinity bfc bfc_trinity \
      rcorrector rcorr_trinity
runseecer: seecer seecer_trinity
subsamp_reads:${SAMP}.subsamp_1.fastq ${SAMP}.subsamp_2.fastq

.DELETE_ON_ERROR:

setup:
	mkdir ${DIR}/error_profiles
	mkdir ${DIR}/assemblies

scripts:
	@echo Downloading Scripts
	mkdir -p ${DIR}/scripts
	cd ${DIR}/scripts && \
	wget https://raw.githubusercontent.com/macmanes/read_error_corr/master/barcodes.fa

download_reads:
	mkdir -p ${DIR}/reads
	cd ${DIR}/reads && \
	curl -LO https://s3.amazonaws.com/macmanes_general/SRR797058_2.fastq.gz && \
	curl -LO https://s3.amazonaws.com/macmanes_general/SRR797058_1.fastq.gz

reference:
	mkdir -p ${DIR}/genome
	cd ${DIR}/genome && \
	wget ftp://ftp.ensembl.org/pub/release-79/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.1.fa.gz && \
	bwa index -p mus Mus_musculus.GRCm38.dna.chromosome.1.fa.gz

${SAMP}.subsamp_1.fastq ${SAMP}.subsamp_2.fastq:
	cd ${DIR}/reads && \
	seqtk sample -s102340 ${READ1} ${SAMP}000000 | sed 's_ H_-H_g' > tee ${SAMP}.subsamp_1.fastq && \
	seqtk sample -s102340 ${READ2} ${SAMP}000000 | sed 's_ H_-H_g' > tee ${SAMP}.subsamp_2.fastq

raw:${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq
	mkdir -p ${DIR}/raw
	cd ${DIR}/raw && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq > ${SAMP}M.raw.sam
	k8 ${BFCDIR}/errstat.js ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.raw.out && \
	mv ${SAMP}M.raw.out ${DIR}/error_profiles/


lighter:${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq
	mkdir -p ${DIR}/lighter${SAMP}M
	cd ${DIR}/lighter${SAMP}M && \
	lighter -K 31 60000000 -r ${DIR}/reads/${SAMP}.subsamp_1.fastq -r ${DIR}/reads/${SAMP}.subsamp_2.fastq -t $(CPU) && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${SAMP}.subsamp_1.cor.fq ${SAMP}.subsamp_2.cor.fq > ${SAMP}M.lighter.sam && \
	k8 ${BFCDIR}/errstat.js ${DIR}/lighter${SAMP}M/${SAMP}M.lighter.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.lighter.out && \
	mv ${SAMP}M.lighter.out ${DIR}/error_profiles/ && \
	rm *sam 

lighter_trinity:
	cd ${DIR}/lighter${SAMP}M && \
	Trinity --seqType fq --max_memory 50G --output ${SAMP}M.trinity.lighter --trimmomatic --left ${SAMP}.subsamp_1.cor.fq --right ${SAMP}.subsamp_2.cor.fq --CPU $(CPU) --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	mv *fasta ${DIR}/assemblies/

bless:
	mkdir -p ${DIR}/bless${SAMP}M
	cd ${DIR}/bless${SAMP}M && \
	mkdir -p ${DIR}/bless${SAMP}M/kmc/bin/ && \
	cp ${KMC}/bin/kmc ${DIR}/bless${SAMP}M/kmc/bin/kmc && \
	mpirun -np $(CPU) bless -notrim -read1 ${DIR}/reads/${SAMP}.subsamp_1.fastq -read2 ${DIR}/reads/${SAMP}.subsamp_2.fastq -prefix ${SAMP}M_bless55 -kmerlength 55 && \
	mpirun -np $(CPU) bless -notrim -read1 ${DIR}/reads/${SAMP}.subsamp_1.fastq -read2 ${DIR}/reads/${SAMP}.subsamp_2.fastq -prefix ${SAMP}M_bless31 -kmerlength 31 && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${SAMP}M_bless55.1.corrected.fastq ${SAMP}M_bless55.2.corrected.fastq > ${SAMP}M.bless55.sam && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${SAMP}M_bless31.1.corrected.fastq ${SAMP}M_bless31.2.corrected.fastq > ${SAMP}M.bless31.sam && \
	k8 ${BFCDIR}/errstat.js ${DIR}/bless${SAMP}M/${SAMP}M.bless55.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bless55.out && \
	k8 ${BFCDIR}/errstat.js ${DIR}/bless${SAMP}M/${SAMP}M.bless31.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bless31.out && \
	mv ${SAMP}M.bless55.out ${DIR}/error_profiles/ && \
	mv ${SAMP}M.bless31.out ${DIR}/error_profiles/ && \
	rm *sam

bless_trinity:
	cd ${DIR}/bless${SAMP}M && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left ${SAMP}M_bless55.1.corrected.fastq --right ${SAMP}M_bless55.2.corrected.fastq --CPU $(CPU) --output ${SAMP}M.trinity_bless55 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left ${SAMP}M_bless31.1.corrected.fastq --right ${SAMP}M_bless31.2.corrected.fastq --CPU $(CPU) --output ${SAMP}M.trinity_bless31 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	mv *fasta ${DIR}/assemblies/

sga:
	mkdir -p ${DIR}/sga${SAMP}M
	cd ${DIR}/sga${SAMP}M && \
	sga preprocess -p 1 ${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq | gzip -1 > out.pe.fq.gz && \
	sga index -a ropebwt -t 8 --no-reverse out.pe.fq.gz && \
	sga correct -t $(CPU) -k 55 --learn out.pe.fq.gz -o sga.55.fq && \
	bwa mem -p -t $(CPU) ${DIR}/genome/mus sga.55.fq > ${SAMP}M.sga55.sam && \
	sga correct -t $(CPU) -k 31 --learn out.pe.fq.gz -o sga.31.fq && \
	bwa mem -p -t $(CPU) ${DIR}/genome/mus sga.31.fq > ${SAMP}M.sga31.sam && \
	k8 ${BFCDIR}/errstat.js ${DIR}/sga${SAMP}M/${SAMP}M.sga55.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.sga55.out && \
	k8 ${BFCDIR}/errstat.js ${DIR}/sga${SAMP}M/${SAMP}M.sga31.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.sga31.out && \
	mv ${SAMP}M.sga31.out ${DIR}/error_profiles/ && \
	mv ${SAMP}M.sga55.out ${DIR}/error_profiles/ && \
	rm *sam

sga_trinity:
	cd ${DIR}/sga${SAMP}M && \
	split-paired-reads.py sga.55.fq && \
	split-paired-reads.py sga.31.fq && \
	Trinity --seqType fq --max_memory 50G --output ${SAMP}M.trinity_sga55 --trimmomatic --left sga.55.fq.1 --right sga.55.fq.2 --CPU $(CPU) --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --seqType fq --max_memory 50G --output ${SAMP}M.trinity_sga31 --trimmomatic --left sga.31.fq.1 --right sga.31.fq.2 --CPU $(CPU) --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	mv *fasta ${DIR}/assemblies/


bfc:
	mkdir -p ${DIR}/bfc${SAMP}M
	cd ${DIR}/bfc${SAMP}M && \
	seqtk mergepe ${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq > inter.fq && \
	bfc -s 50m -k55 -t $(CPU) inter.fq | tee bfc55.corr.fq | bwa mem -p -t $(CPU) ${DIR}/genome/mus - > ${SAMP}M.bfc55.sam && \
	bfc -s 50m -k31 -t $(CPU) inter.fq | tee bfc31.corr.fq | bwa mem -p -t $(CPU) ${DIR}/genome/mus - > ${SAMP}M.bfc31.sam && \
	k8 ${BFCDIR}/errstat.js ${DIR}/bfc${SAMP}M/${SAMP}M.bfc55.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bfc55.out && \
	k8 ${BFCDIR}/errstat.js ${DIR}/bfc${SAMP}M/${SAMP}M.bfc31.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bfc31.out && \
	mv ${SAMP}M.bfc31.out ${DIR}/error_profiles/ && \
	mv ${SAMP}M.bfc55.out ${DIR}/error_profiles/ && \
	rm *sam

bfc_trinity:
	cd ${DIR}/bfc${SAMP}M && \
	split-paired-reads.py bfc55.corr.fq && \
	split-paired-reads.py bfc31.corr.fq && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left bfc55.corr.fq.1 --right bfc55.corr.fq.2 --CPU $(CPU) --output ${SAMP}M.trinity_bfc55 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left bfc31.corr.fq.1 --right bfc31.corr.fq.2 --CPU $(CPU) --output ${SAMP}M.trinity_bfc31 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	mv *fasta ${DIR}/assemblies/


seecer:
	mkdir -p ${DIR}/seecer${SAMP}M
	cd ${DIR}/seecer${SAMP}M && \
	run_seecer.sh -t . -k 31 ${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${DIR}/reads/${SAMP}.subsamp_1.fastq_corrected.fa ${DIR}/reads/${SAMP}.subsamp_2.fastq_corrected.fa > ${SAMP}M.seecer.sam && \
	k8 ${BFCDIR}/errstat.js ${DIR}/seecer${SAMP}M/${SAMP}M.seecer.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.seecer.out && \
	mv ${SAMP}M.seecer.out ${DIR}/error_profiles/ && \
	rm *sam


seecer_trinity:
	cd ${DIR}/seecer${SAMP}M && \
	Trinity --seqType fa --max_memory 50G --output ${SAMP}M.seecer --trimmomatic --left ${DIR}/reads/${SAMP}.subsamp_1.fastq_corrected.fa --right ${DIR}/reads/${SAMP}.subsamp_2.fastq_corrected.fa --CPU $(CPU) --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	mv *fasta ${DIR}/assemblies/

rcorrector:
	mkdir -p ${DIR}/rcorr${SAMP}M
	cd ${DIR}/rcorr${SAMP}M && \
	perl ~/Rcorrector/run_rcorrector.pl -t $(CPU) -k 25 -1 ${DIR}/reads/${SAMP}.subsamp_1.fastq -2 ${DIR}/reads/${SAMP}.subsamp_2.fastq && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${DIR}/rcorr${SAMP}M/${SAMP}.subsamp_1.cor.fq ${DIR}/rcorr${SAMP}M/${SAMP}.subsamp_2.cor.fq > ${SAMP}M.rcorr.sam && \
	k8 ${BFCDIR}/errstat.js ${DIR}/rcorr${SAMP}M/${SAMP}M.rcorr.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.rcorr.out && \
	mv ${SAMP}M.rcorr.out ${DIR}/error_profiles/ && \
	rm *sam

rcorr_trinity:
	cd ${DIR}/rcorr${SAMP}M && \
	Trinity --seqType fq --output ${SAMP}M.trinity_rcorr --max_memory 50G --trimmomatic --left ${DIR}/rcorr${SAMP}M/${SAMP}.subsamp_1.cor.fq --right ${DIR}/rcorr${SAMP}M/${SAMP}.subsamp_2.cor.fq --CPU $(CPU) --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	mv *fasta ${DIR}/assemblies/


trinity_raw:${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq
	mkdir -p ${DIR}/trinity_${SAMP}M
	cd ${DIR}/trinity_${SAMP}M && \
	Trinity --seqType fq --max_memory 50G --trimmomatic --left $< --right $(word 2,$^) --CPU $(CPU) --output trinity_${SAMP}M.P2.raw --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	mv *fasta ${DIR}/assemblies/
