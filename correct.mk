#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#	./correct.mk all CPU=16 SAMP=10
#
#
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
SAMP=10
#the number of reads in millions to subsample.


all: scripts download_reads subsamp_reads reference raw lighter bless sga bfc seecer rcorrector stats trinity_bfc trinity_raw


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
subsamp_reads:
	cd ${DIR}/reads && \
	seqtk sample -s102340 SRR797058_1.fastq.gz ${SAMP}000000 > subsamp_1.fastq && \
	seqtk sample -s102340 SRR797058_2.fastq.gz ${SAMP}000000 > subsamp_2.fastq && \
	sed -i 's_ H_-H_g' subsamp_{1,2}.fastq

reference:
	mkdir -p ${DIR}/genome
	cd ${DIR}/genome && \
	wget ftp://ftp.ensembl.org/pub/release-79/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.1.fa.gz && \
	bwa index -p mus Mus_musculus.GRCm38.dna.chromosome.1.fa.gz

raw:${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq
	mkdir -p ${DIR}/raw
	cd ${DIR}/raw && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq > ${SAMP}M.raw.sam && \
	Trinity --seqType fq --max_memory 10G --trimmomatic --left $< --right $(word 2,$^) --CPU $(CPU) --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"


lighter:${DIR}/lighter${SAMP}M/subsamp_1.cor.fq ${DIR}/lighter${SAMP}M/subsamp_2.cor.fq
	mkdir -p ${DIR}/lighter${SAMP}M
	cd ${DIR}/lighter${SAMP}M && \
	lighter -K 31 60000000 -r ${DIR}/reads/subsamp_1.fastq -r ${DIR}/reads/subsamp_2.fastq -t $(CPU) && \
	bwa mem -t $(CPU) ${DIR}/genome/mus subsamp_1.cor.fq subsamp_2.cor.fq > ${SAMP}M.lighter.sam && \
	Trinity --seqType fq --max_memory 10G --trimmomatic --left $< --right $(word 2,$^) --CPU $(CPU) --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"


bless:
	mkdir -p ${DIR}/bless${SAMP}M
	mkdir -p ${DIR}/bless${SAMP}M/kmc/bin/
	cp /home/ubuntu/v0p24/kmc/bin/kmc /mnt/bless${SAMP}M/kmc/bin/kmc
	cd ${DIR}/bless${SAMP}M && \
	mpirun -np $(CPU) bless -read1 ${DIR}/reads/subsamp_1.fastq -read2 ${DIR}/reads/subsamp_2.fastq -prefix ${SAMP}M_bless55 -kmerlength 55 && \
	mpirun -np $(CPU) bless -read1 ${DIR}/reads/subsamp_1.fastq -read2 ${DIR}/reads/subsamp_2.fastq -prefix ${SAMP}M_bless33 -kmerlength 33 && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${SAMP}M_bless55.1.corrected.fastq ${SAMP}M_bless55.2.corrected.fastq > ${SAMP}M.bless55.sam && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${SAMP}M_bless33.1.corrected.fastq ${SAMP}M_bless33.2.corrected.fastq > ${SAMP}M.bless33.sam


sga:
	mkdir -p ${DIR}/sga${SAMP}M
	cd ${DIR}/sga${SAMP}M && \
	sga preprocess -p 1 ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq | gzip -1 > out.pe.fq.gz && \
	sga index -a ropebwt -t 8 --no-reverse out.pe.fq.gz && \
	sga correct -t $(CPU) -k 55 --learn out.pe.fq.gz -o sga.55.fq && \
	sga correct -t $(CPU) -k 33 --learn out.pe.fq.gz -o sga.33.fq && \
	bwa mem -p -t $(CPU) ${DIR}/genome/mus sga.55.fq > ${SAMP}M.sga55.sam && \
	bwa mem -p -t $(CPU) ${DIR}/genome/mus sga.33.fq > ${SAMP}M.sga33.sam

bfc:
	mkdir -p ${DIR}/bfc${SAMP}M
	cd ${DIR}/bfc${SAMP}M && \
	seqtk mergepe ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq > inter.fq && \
	bfc -s 50m -k55 -t $(CPU) inter.fq > bfc55.corr.fq && \
	bfc -s 50m -k33 -t $(CPU) inter.fq > bfc33.corr.fq && \
	bwa mem -p -t $(CPU) ${DIR}/genome/mus bfc55.corr.fq > ${SAMP}M.bfc55.sam && \
	bwa mem -p -t $(CPU) ${DIR}/genome/mus bfc33.corr.fq > ${SAMP}M.bfc33.sam

seecer:${DIR}/reads/subsamp_1.fastq_corrected.fa ${DIR}/reads/subsamp_2.fastq_corrected.fa
	mkdir -p ${DIR}/seecer${SAMP}M
	cd ${DIR}/seecer${SAMP}M && \
	run_seecer.sh -t . -k 31 ${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq && \
	bwa mem -t $(CPU) ${DIR}/genome/mus ${DIR}/reads/subsamp_1.fastq_corrected.fa ${DIR}/reads/subsamp_2.fastq_corrected.fa > ${SAMP}M.seecer.sam && \
	Trinity --seqType fa --max_memory 10G --trimmomatic --left $< --right $(word 2,$^) --CPU $(CPU) --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"


rcorrector:
	mkdir -p ${DIR}/rcorr${SAMP}M
	cd ${DIR}/rcorr${SAMP}M && \
	perl /share/Rcorrector/run_rcorrector.pl -t $(CPU) -k 25 -1 ${DIR}/reads/subsamp_1.fastq -2 ${DIR}/reads/subsamp_2.fastq
stats:
	mkdir -p ${DIR}/stats${SAMP}M
	cd ${DIR}/stats${SAMP}M && \
	k8 ~/bfc/errstat.js ${DIR}/bfc${SAMP}M/${SAMP}M.bfc55.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bfc55.out && \
	k8 ~/bfc/errstat.js ${DIR}/bfc${SAMP}M/${SAMP}M.bfc33.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bfc33.out && \
	k8 ~/bfc/errstat.js ${DIR}/sga${SAMP}M/${SAMP}M.sga55.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.sga33.out && \
	k8 ~/bfc/errstat.js ${DIR}/sga${SAMP}M/${SAMP}M.sga33.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.sga33.out && \
	k8 ~/bfc/errstat.js ${DIR}/bless${SAMP}M/${SAMP}M.bless33.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bless33.out && \
	k8 ~/bfc/errstat.js ${DIR}/bless${SAMP}M/${SAMP}M.bless33.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bless33.out && \
	k8 ~/bfc/errstat.js ${DIR}/lighter${SAMP}M/${SAMP}M.lighter.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.lighter.out && \
	k8 ~/bfc/errstat.js ${DIR}/seecer${SAMP}M/${SAMP}M.seecer.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.seecer.out


trinity_raw:${DIR}/reads/subsamp_1.fastq ${DIR}/reads/subsamp_2.fastq
	mkdir -p ${DIR}/trinity_${SAMP}M
	cd ${DIR}/trinity_${SAMP}M && \
	Trinity --seqType fq --max_memory 10G --trimmomatic \
	--left $< \
	--right $(word 2,$^) \
	--CPU $(CPU) --output trinity_${SAMP}M.P2.raw --inchworm_cpu 10 --full_cleanup \
	--quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25"
