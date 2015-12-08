#!/usr/bin/make -rRsf

SHELL=/bin/bash -o pipefail

#USAGE:
#
#	for i in 1 2 5 10 20 40 60 80 100; do ./correct.mk main SAMP=$i CPU=36; done
#

MAKEDIR := $(dir $(firstword $(MAKEFILE_LIST)))
DIR := ${CURDIR}
CPU=16
SAMP=10
BFC ?= ${shell which bfc}
BFCDIR := $(dir $(firstword $(BFC)))
BLESS ?= ${shell which bless}
KMC := $(dir $(firstword $(BLESS)))
RCORR ?= ${shell which rcorrector}
RCORRDIR := $(dir $(firstword $(RCORR)))
READ1=SRR1659968_1.fastq
READ2=SRR1659968_2.fastq


prep: download_reads scripts reference
main: setup subsamp_reads raw trinity_raw lighter lighter_trinity bless bless_trinity bfc bfc_trinity rcorrector rcorr_trinity seecer seecer_trinity
subsamp_reads:${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq

.DELETE_ON_ERROR:

setup:
	mkdir -p ${DIR}/error_profiles
	mkdir -p ${DIR}/assemblies

scripts:
	@echo Downloading Scripts
	mkdir -p ${DIR}/scripts
	cd ${DIR}/scripts && \
	wget https://raw.githubusercontent.com/macmanes/read_error_corr/master/barcodes.fa

download_reads:
	mkdir -p ${DIR}/reads
	cd ${DIR}/reads && \
	curl -LO curl -LO ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR165/008/SRR1659968/SRR1659968_1.fastq.gz && \
	curl -LO curl -LO ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR165/008/SRR1659968/SRR1659968_2.fastq.gz && \
	gzip -d *gz

reference:
	mkdir -p ${DIR}/genome
	cd ${DIR}/genome && \
	curl -LO http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz && \
	tar -zxf hg38.chromFa.tar.gz && \
	bwa index -p homo chroms/chr1.fa

${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq:
	cd ${DIR}/reads && \
	seqtk sample -s102340 ${READ1} ${SAMP}000000 | sed 's_ _-_g' > ${SAMP}.subsamp_1.fastq && \
	seqtk sample -s102340 ${READ2} ${SAMP}000000 | sed 's_ _-_g' > ${SAMP}.subsamp_2.fastq

raw:
	mkdir -p ${DIR}/raw
	cd ${DIR}/raw && \
	bwa mem -t $(CPU) ${DIR}/genome/homo ${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq > ${SAMP}M.raw.sam
	k8 ${BFCDIR}/errstat.js ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.raw.out && \
	mv ${SAMP}M.raw.out ${DIR}/error_profiles/


lighter:${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq
	mkdir -p ${DIR}/lighter${SAMP}M
	cd ${DIR}/lighter${SAMP}M && \
	lighter -K 31 60000000 -r ${DIR}/reads/${SAMP}.subsamp_1.fastq -r ${DIR}/reads/${SAMP}.subsamp_2.fastq -t $(CPU) && \
	bwa mem -t $(CPU) ${DIR}/genome/homo ${SAMP}.subsamp_1.cor.fq ${SAMP}.subsamp_2.cor.fq > ${SAMP}M.lighter.sam && \
	k8 ${BFCDIR}/errstat.js ${DIR}/lighter${SAMP}M/${SAMP}M.lighter.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.lighter.out && \
	mv ${SAMP}M.lighter.out ${DIR}/error_profiles/ && \
	rm *sam

lighter_trinity:
	cd ${DIR}/lighter${SAMP}M && \
	Trinity --SS_lib_type RF --seqType fq --max_memory 50G --output ${SAMP}M.trinity.lighter --trimmomatic --left ${SAMP}.subsamp_1.cor.fq --right ${SAMP}.subsamp_2.cor.fq --CPU $(CPU) --inchworm_cpu 10 --full_cleanup --quality_trimmi
ng_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o ${SAMP}M.lighter -g ${SAMP}M.trinity.lighter.Trinity.fasta -m Trans --cpu $(CPU) -l ~/BUSCO_v1.1b1/vertebrata && \
	transrate -o ${SAMP}M.lighter -a ${SAMP}M.trinity.lighter.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	mv *fasta ${DIR}/assemblies/

bless:
	mkdir -p ${DIR}/bless${SAMP}M
	cd ${DIR}/bless${SAMP}M && \
	mkdir -p ${DIR}/bless${SAMP}M/kmc/bin/ && \
	cp ${KMC}/kmc/bin/kmc ${DIR}/bless${SAMP}M/kmc/bin/kmc && \
	mpirun -np $(CPU) bless -notrim -read1 ${DIR}/reads/${SAMP}.subsamp_1.fastq -read2 ${DIR}/reads/${SAMP}.subsamp_2.fastq -prefix ${SAMP}M_bless55 -kmerlength 55 && \
	mpirun -np $(CPU) bless -notrim -read1 ${DIR}/reads/${SAMP}.subsamp_1.fastq -read2 ${DIR}/reads/${SAMP}.subsamp_2.fastq -prefix ${SAMP}M_bless31 -kmerlength 31 && \
	bwa mem -t $(CPU) ${DIR}/genome/homo ${SAMP}M_bless55.1.corrected.fastq ${SAMP}M_bless55.2.corrected.fastq > ${SAMP}M.bless55.sam && \
	bwa mem -t $(CPU) ${DIR}/genome/homo ${SAMP}M_bless31.1.corrected.fastq ${SAMP}M_bless31.2.corrected.fastq > ${SAMP}M.bless31.sam && \
	k8 ${BFCDIR}/errstat.js ${DIR}/bless${SAMP}M/${SAMP}M.bless55.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bless55.out && \
	k8 ${BFCDIR}/errstat.js ${DIR}/bless${SAMP}M/${SAMP}M.bless31.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bless31.out && \
	mv ${SAMP}M.bless55.out ${DIR}/error_profiles/ && \
	mv ${SAMP}M.bless31.out ${DIR}/error_profiles/ && \
	rm *sam

bless_trinity:
	cd ${DIR}/bless${SAMP}M && \
	Trinity --SS_lib_type RF --seqType fq --max_memory 50G --trimmomatic --left ${SAMP}M_bless55.1.corrected.fastq --right ${SAMP}M_bless55.2.corrected.fastq --CPU $(CPU) --output ${SAMP}M.trinity_bless55 --inchworm_cpu 10 --full_clea
nup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --SS_lib_type RF --seqType fq --max_memory 50G --trimmomatic --left ${SAMP}M_bless31.1.corrected.fastq --right ${SAMP}M_bless31.2.corrected.fastq --CPU $(CPU) --output ${SAMP}M.trinity_bless31 --inchworm_cpu 10 --full_clea
nup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o ${SAMP}M.bless55  -g ${SAMP}M.trinity_bless55.Trinity.fasta -m Trans --cpu $(CPU) -l ~/BUSCO_v1.1b1/vertebrata && \
	python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o ${SAMP}M.bless31 -g ${SAMP}M.trinity_bless31.Trinity.fasta -m Trans --cpu $(CPU) -l ~/BUSCO_v1.1b1/vertebrata && \
	transrate -o ${SAMP}M.bless55 -a ${SAMP}M.trinity_bless55.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	transrate -o ${SAMP}M.bless31 -a ${SAMP}M.trinity_bless31.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	mv *fasta ${DIR}/assemblies/

bfc ${DIR}/reads/${SAMP}.inter.fq:
	mkdir -p ${DIR}/bfc${SAMP}M
	cd ${DIR}/bfc${SAMP}M && \
	seqtk mergepe ${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq > ${DIR}/reads/${SAMP}.inter.fq && \
	bfc -s 50m -k55 -t $(CPU) ${DIR}/reads/${SAMP}.inter.fq | tee bfc55.corr.fq | bwa mem -p -t $(CPU) ${DIR}/genome/homo - > ${SAMP}M.bfc55.sam && \
	bfc -s 50m -k31 -t $(CPU) ${DIR}/reads/${SAMP}.inter.fq | tee bfc31.corr.fq | bwa mem -p -t $(CPU) ${DIR}/genome/homo - > ${SAMP}M.bfc31.sam && \
	k8 ${BFCDIR}/errstat.js ${DIR}/bfc${SAMP}M/${SAMP}M.bfc55.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bfc55.out && \
	k8 ${BFCDIR}/errstat.js ${DIR}/bfc${SAMP}M/${SAMP}M.bfc31.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.bfc31.out && \
	mv ${SAMP}M.bfc31.out ${DIR}/error_profiles/ && \
	mv ${SAMP}M.bfc55.out ${DIR}/error_profiles/ && \
	rm *sam

bfc_trinity:
	cd ${DIR}/bfc${SAMP}M && \
	split-paired-reads.py bfc55.corr.fq && \
	split-paired-reads.py bfc31.corr.fq && \
	Trinity --SS_lib_type RF --seqType fq --max_memory 50G --trimmomatic --left bfc55.corr.fq.1 --right bfc55.corr.fq.2 --CPU $(CPU) --output ${SAMP}M.trinity_bfc55 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINAC
LIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --SS_lib_type RF --seqType fq --max_memory 50G --trimmomatic --left bfc31.corr.fq.1 --right bfc31.corr.fq.2 --CPU $(CPU) --output ${SAMP}M.trinity_bfc31 --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINAC
LIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o ${SAMP}M.bfc55 -g ${SAMP}M.trinity_bfc55.Trinity.fasta -m Trans --cpu $(CPU) -l ~/BUSCO_v1.1b1/vertebrata && \
	python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o ${SAMP}M.bfc31 -g ${SAMP}M.trinity_bfc31.Trinity.fasta -m Trans --cpu $(CPU) -l ~/BUSCO_v1.1b1/vertebrata && \
	transrate -o ${SAMP}M.bfc55 -a ${SAMP}M.trinity_bfc55.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	transrate -o ${SAMP}M.bfc31 -a ${SAMP}M.trinity_bfc31.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	mv *fasta ${DIR}/assemblies/


seecer:
	mkdir -p ${DIR}/seecer${SAMP}M
	cd ${DIR}/seecer${SAMP}M && \
	run_seecer.sh -t . -k 31 ${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq && \
	bwa mem -t $(CPU) ${DIR}/genome/homo ${DIR}/reads/${SAMP}.subsamp_1.fastq_corrected.fa ${DIR}/reads/${SAMP}.subsamp_2.fastq_corrected.fa > ${SAMP}M.seecer.sam && \
	k8 ${BFCDIR}/errstat.js ${DIR}/seecer${SAMP}M/${SAMP}M.seecer.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.seecer.out && \
	mv ${SAMP}M.seecer.out ${DIR}/error_profiles/ && \
	rm *sam


seecer_trinity:
	cd ${DIR}/seecer${SAMP}M && \
	Trinity --seqType fa --SS_lib_type RF --max_memory 50G --output ${SAMP}M.seecer --trimmomatic --left ${DIR}/reads/${SAMP}.subsamp_1.fastq_corrected.fa --right ${DIR}/reads/${SAMP}.subsamp_2.fastq_corrected.fa --CPU $(CPU) --inchwo
rm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o ${SAMP}M.seecer -g ${SAMP}M.trinity_seecer.Trinity.fasta -m Trans --cpu $(CPU) -l ~/BUSCO_v1.1b1/vertebrata && \
	transrate -o ${SAMP}M.seecer -a ${SAMP}M.trinity_seecer.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	mv *fasta ${DIR}/assemblies/

rcorrector:${DIR}/reads/${SAMP}.inter.fq
	mkdir -p ${DIR}/rcorr${SAMP}M
	cd ${DIR}/rcorr${SAMP}M && \
	perl ${RCORRDIR}/run_rcorrector.pl -t $(CPU) -k 31 -i ${DIR}/reads/${SAMP}.inter.fq -stdout | sed 1,8d | tee rcorr31.corr.fq| bwa mem -t $(CPU) ${DIR}/genome/homo - > ${SAMP}M.rcorr31.sam && \
	perl ${RCORRDIR}/run_rcorrector.pl -t $(CPU) -k 55 -i ${DIR}/reads/${SAMP}.inter.fq -stdout | sed 1,8d | tee rcorr55.corr.fq | bwa mem -t $(CPU) ${DIR}/genome/homo - > ${SAMP}M.rcorr55.sam && \
	k8 ${BFCDIR}/errstat.js ${SAMP}M.rcorr31.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.rcorr31.out && \
	k8 ${BFCDIR}/errstat.js ${SAMP}M.rcorr55.sam ${DIR}/raw/${SAMP}M.raw.sam | tail -11 > ${SAMP}M.rcorr55.out && \
	mv ${SAMP}M.*.out ${DIR}/error_profiles/ && \
	rm *sam

rcorr_trinity:
	cd ${DIR}/rcorr${SAMP}M && \
	split-paired-reads.py rcorr31.corr.fq && \
	split-paired-reads.py rcorr55.corr.fq && \
	Trinity --seqType fq --SS_lib_type RF --output ${SAMP}M.trinity_rcorr31 --max_memory 50G --trimmomatic --left ${DIR}/rcorr${SAMP}M/rcorr31.corr.fq.1 --right ${DIR}/rcorr${SAMP}M/rcorr31.corr.fq.2 --CPU $(CPU) --inchworm_cpu 10 --f
ull_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	Trinity --seqType fq --SS_lib_type RF --output ${SAMP}M.trinity_rcorr55 --max_memory 50G --trimmomatic --left ${DIR}/rcorr${SAMP}M/rcorr55.corr.fq.1 --right ${DIR}/rcorr${SAMP}M/rcorr55.corr.fq.2 --CPU $(CPU) --inchworm_cpu 10 --f
ull_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scripts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o ${SAMP}M.rcorr55 -g ${SAMP}M.trinity_rcorr55.Trinity.fasta -m Trans --cpu $(CPU) -l ~/BUSCO_v1.1b1/vertebrata && \
	python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o ${SAMP}M.rcocc31 -g ${SAMP}M.trinity_rcorr31.Trinity.fasta -m Trans --cpu $(CPU) -l ~/BUSCO_v1.1b1/vertebrata && \
	transrate -o ${SAMP}M.rcorr55 -a ${SAMP}M.trinity_rcorr55.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	transrate -o ${SAMP}M.rcorr31 -a ${SAMP}M.trinity_rcorr31.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	mv *fasta ${DIR}/assemblies/


trinity_raw:${DIR}/reads/${SAMP}.subsamp_1.fastq ${DIR}/reads/${SAMP}.subsamp_2.fastq
	mkdir -p ${DIR}/trinity_${SAMP}M
	cd ${DIR}/trinity_${SAMP}M && \
	Trinity --seqType fq --max_memory 50G --SS_lib_type RF --trimmomatic --left $< --right $(word 2,$^) --CPU $(CPU) --output trinity_${SAMP}M.P2.raw --inchworm_cpu 10 --full_cleanup --quality_trimming_params "ILLUMINACLIP:${DIR}/scri
pts/barcodes.fa:2:40:15 LEADING:2 TRAILING:2 MINLEN:25" && \
	python3 ~/BUSCO_v1.1b1/BUSCO_v1.1b1.py -o ${SAMP}M.raw -g trinity_${SAMP}M.P2.raw.Trinity.fasta -m Trans --cpu $(CPU) -l ~/BUSCO_v1.1b1/vertebrata && \
	transrate -o ${SAMP}M.raw -a trinity_${SAMP}M.P2.raw.Trinity.fasta --left ${DIR}/reads/${READ1} --right ${DIR}/reads/${READ2} -t $(CPU) && \
	mv *fasta ${DIR}/assemblies/
