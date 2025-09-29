#!/usr/bin/env bash

# set -x
# set -v
# set -u
set -e
set -o pipefail

# conda activation hook
eval "$(conda shell.bash hook)"

# set common variables
BASEDIR=/home/bioinf/Projects/mutants/2023-03-09

# DEF SMARTDENOVO DRAFT ASSEMBLE
smartdenovo_function() {

    # set variables
    local CURRENT_SUBSAMPLE=$1   # 'local' makes the variable visible only inside the function
    
    smartdenovo.pl -t $THREADS -p sample${CURRENT_SUBSAMPLE}.tmd.flt -c 1 -J $MIN_READ_LENGTH $BASEDIR/$CURRENT_DIR/longreads/subsamples/sample_${CURRENT_SUBSAMPLE}.fastq > sample${CURRENT_SUBSAMPLE}.tmd.flt.mak
    make -f sample${CURRENT_SUBSAMPLE}.tmd.flt.mak
    mv sample${CURRENT_SUBSAMPLE}.tmd.flt.dmo.cns $BASEDIR/$CURRENT_DIR/assemblies/sample${CURRENT_SUBSAMPLE}.smdnv.fasta
    # the folder can be deleted   
}

# DEF CANU DRAFT ASSEMBLE
canu_function() {
    
    # set variables
    local CURRENT_SUBSAMPLE=$1
    
    canu -d sample${CURRENT_SUBSAMPLE} -p sample${CURRENT_SUBSAMPLE}.tmd.flt genomeSize=30000 minReadLength=$MIN_READ_LENGTH  minOverlapLength=100 -nanopore ${BASEDIR}/${CURRENT_DIR}/longreads/subsamples/sample_${CURRENT_SUBSAMPLE}.fastq
    # minOverlapLength value must not be greater than minReadLength value (default 500)
    mv sample${CURRENT_SUBSAMPLE}/sample${CURRENT_SUBSAMPLE}.tmd.flt.contigs.fasta ${BASEDIR}/${CURRENT_DIR}/assemblies/sample${CURRENT_SUBSAMPLE}.canu.fasta
    # the folder can be deleted   
}

# DEF NECAT DRAFT ASSEMBLE
necat_function() {

    # set variables
    local CURRENT_SUBSAMPLE=$1
    
    sed -i 's/sample.*/sample'$CURRENT_SUBSAMPLE'/g' /home/bioinf/necat/template_config.txt
    sed -i 's/MIN_READ_LENGTH=.*/MIN_READ_LENGTH='$MIN_READ_LENGTH'/g' /home/bioinf/necat/template_config.txt
    sed -i 's/THREADS=.*/THREADS='$THREADS'/g' /home/bioinf/necat/template_config.txt
    echo "$BASEDIR/$CURRENT_DIR/longreads/subsamples/sample_${CURRENT_SUBSAMPLE}.fastq" >| /home/bioinf/necat/reads_list.txt
    necat correct /home/bioinf/necat/template_config.txt
    necat assemble /home/bioinf/necat/template_config.txt
    mv sample${CURRENT_SUBSAMPLE}/4-fsa/polished_contigs.fasta $BASEDIR/$CURRENT_DIR/assemblies/sample${CURRENT_SUBSAMPLE}.necat.fasta
    # the folder can be deleted   
}

# DEF FLYE DRAFT ASSEMBLE
flye_function() {    

    # set variables
    local CURRENT_SUBSAMPLE=$1
    
    conda activate flye
    flye --genome-size 30000 --threads $THREADS --nano-hq ${BASEDIR}/${CURRENT_DIR}/longreads/subsamples/sample_${CURRENT_SUBSAMPLE}.fastq --out-dir sample${CURRENT_SUBSAMPLE}
    mv sample${CURRENT_SUBSAMPLE}/assembly.fasta $BASEDIR/$CURRENT_DIR/assemblies/sample${CURRENT_SUBSAMPLE}.flye.fasta
    conda deactivate    
}

# DEF RECONCILING
reconciling() {

    #set variables
    local CURRENT_BARCODE=0$1
    local THREADS=$2
    local CURRENT_DIR=barcode$CURRENT_BARCODE
    local CURRENT_SUFFIX=bc$CURRENT_BARCODE
    
    # TRYCYCLER CLUSTER RECONCILING
    conda activate progs
    cd $BASEDIR/$CURRENT_DIR
    trycycler reconcile --threads $THREADS --reads longreads/${CURRENT_SUFFIX}.tmd.flt.fastq.gz --cluster_dir trycycler/cluster_001
    
    # TRYCYCLER CREATING CLUSTER MSA
    trycycler msa --threads $THREADS --cluster_dir trycycler/cluster_001
    
    # TRYCYCLER CLUSTER PARTITIONING
    trycycler partition --threads $THREADS --reads longreads/${CURRENT_SUFFIX}.tmd.flt.fastq.gz --cluster_dirs trycycler/cluster_*
    
    # TRYCYCLER CONCENSUS GENERATING
    trycycler consensus --linear --threads $THREADS --cluster_dir trycycler/cluster_001
    
    # HOMOPOLISH CONSENSUS POLISHING
    conda activate homopolish
    cd /home/bioinf/homopolish
    python3 homopolish.py polish --threads $THREADS --sketch_path virus.msh --model_path R9.4.pkl --assembly $BASEDIR/$CURRENT_DIR/trycycler/cluster_001/7_final_consensus.fasta --output_dir $BASEDIR/$CURRENT_DIR/homopolish
    conda deactivate
    
    # MEDAKA CONSENSUS POLISHING
    conda activate medaka
    cd $BASEDIR/$CURRENT_DIR
    medaka_consensus -t $THREADS -m r941_min_hac_g507 -i kraken2/${CURRENT_SUFFIX}_taxid-filtered.tmd.fastq.gz -d homopolish/7_final_consensus_homopolished.fasta -o medaka
    cp medaka/consensus.fasta ${CURRENT_SUFFIX}.consensus.tmd.flt.pol.fasta
    conda deactivate
    
    # FINAL QUALITY CONTROL
    # igv quality control
    conda activate progs
    minimap2 -ax map-ont -t $THREADS ${CURRENT_SUFFIX}.consensus.tmd.flt.pol.fasta kraken2/${CURRENT_SUFFIX}_taxid-filtered.tmd.fastq.gz | samtools sort > ${CURRENT_SUFFIX}.tmd.flt.pol.fasta.bam
    samtools index ${CURRENT_SUFFIX}.tmd.flt.pol.fasta.bam
    
    # polished assembly length
    grep -v ">" ${CURRENT_SUFFIX}.consensus.tmd.flt.pol.fasta | wc | awk '{print "Polished assembly length is " $3-$1}'
    conda deactivate
}

# DEF PIPELINE
pipeline() {

    local CURRENT_BARCODE=$1
    local THREADS=$2
    local MIN_READ_LENGTH=$3 # Flye по дефолту не работает с ридами длиной менее 1000bp
    local CURRENT_DIR=barcode$CURRENT_BARCODE
    local CURRENT_SUFFIX=bc$CURRENT_BARCODE
    
    # PORECHOP TRIMMING
    conda activate porechop
    porechop -i $BASEDIR/$CURRENT_DIR --format fastq.gz -b $BASEDIR/$CURRENT_DIR/porechop
    zcat $BASEDIR/$CURRENT_DIR/porechop/BC${CURRENT_BARCODE}.fastq.gz $BASEDIR/$CURRENT_DIR/porechop/none.fastq.gz | gzip > $BASEDIR/$CURRENT_DIR/porechop/${CURRENT_SUFFIX}.mrg.fastq.gz
    conda deactivate
    
    # KRAKEN2 CONTAMINATION CONTROL
    # kraken2 analysis
    conda activate progs
    mkdir $BASEDIR/$CURRENT_DIR/kraken2
    cd $BASEDIR/$CURRENT_DIR/kraken2
    export KRAKEN2_DB_PATH=~/kraken2 # set kraken2 path
    kraken2 --db standard --threads $THREADS --use-names --report ${CURRENT_SUFFIX}.report $BASEDIR/$CURRENT_DIR/porechop/${CURRENT_SUFFIX}.mrg.fastq.gz > ${CURRENT_SUFFIX}.kraken2.out
    conda deactivate
    
    # manual kraken2 visualization here:
    # sudo docker run --rm -p 5000:80 florianbw/pavian # Pavian launch on http://147.8.70.84:5000/
    
    # taxid-based reads binning
    conda activate krakentools
    read -p "$(tput setaf 3)Please enter taxid: $(tput sgr 0)" TAXID
    extract_kraken_reads.py -k ${CURRENT_SUFFIX}.kraken2.out -s $BASEDIR/$CURRENT_DIR/porechop/${CURRENT_SUFFIX}.mrg.fastq.gz -t $TAXID --include-children -r ${CURRENT_SUFFIX}.report --fastq-output --output ${CURRENT_SUFFIX}_taxid-filtered.tmd.fastq
    gzip ${CURRENT_SUFFIX}_taxid-filtered.tmd.fastq
    conda deactivate

    # FASTQC QUALITY CONTROL
    # fastqc running
    conda activate progs
    mkdir $BASEDIR/$CURRENT_DIR/fastqc  
    fastqc $BASEDIR/$CURRENT_DIR/kraken2/${CURRENT_SUFFIX}_taxid-filtered.tmd.fastq.gz -q -t $THREADS -o $BASEDIR/$CURRENT_DIR/fastqc
    # quality checkpoint
    read -p "$(tput setaf 3)Please check the FASTQC results. Continue? [y/n] $(tput sgr 0)" YN
    echo $YN
    if [ $YN != y ]
    then
        echo "$(tput setaf 3)Stopping pipeline$(tput sgr 0)"
        exit 0
    fi

    # LONG READS FILTERING
    mkdir $BASEDIR/$CURRENT_DIR/longreads
    cd $BASEDIR/$CURRENT_DIR/longreads
    filtlong --min_length $MIN_READ_LENGTH --keep_percent 95 $BASEDIR/$CURRENT_DIR/kraken2/${CURRENT_SUFFIX}_taxid-filtered.tmd.fastq.gz | gzip > ${CURRENT_SUFFIX}.tmd.flt.fastq.gz

    # TRYCYCLER DATA SUBSAMPLING
    trycycler subsample --genome_size 30000 --count 1 --reads ${CURRENT_SUFFIX}.tmd.flt.fastq.gz --out_dir subsamples
    
    # DRAFT ASSEMBLIES
    # running assemblers
    mkdir $BASEDIR/$CURRENT_DIR/assemblies

    # running smartdenovo assembler
    mkdir $BASEDIR/$CURRENT_DIR/smartdenovo
    cd $BASEDIR/$CURRENT_DIR/smartdenovo
    
    for s in $(seq -f "%02g" 1 1)
    do
        smartdenovo_function $s
    done
    
    # running canu assembler
    mkdir /home/bioinf/canu/phages_run1_${CURRENT_SUFFIX}
    cd /home/bioinf/canu/phages_run1_${CURRENT_SUFFIX}
    
    for s in $(seq -f "%02g" 1 1)
    do
        canu_function $s
    done
    
    # running necat assembler
    mkdir /home/bioinf/necat/phages_run1_${CURRENT_SUFFIX}
    cd /home/bioinf/necat/phages_run1_${CURRENT_SUFFIX}
    
    for s in $(seq -f "%02g" 1 5)
    do
        necat_function $s
    done
    
    # running flye assembler
    mkdir $BASEDIR/$CURRENT_DIR/flye
    cd $BASEDIR/$CURRENT_DIR/flye
    
    for s in $(seq -f "%02g" 1 1)
    do
        flye_function $s
    done
    
    # quality checkpoint
    # removing zero-length and trash assemblies (e.g. Bandage can be used here)
    read -p "$(tput setaf 3)Please examine assemblies and remove bad assemblies. Continue? [y/n] $(tput sgr 0)" YN
    echo $YN
    if [ $YN != y ]
    then
        echo "$(tput setaf 3)Stopping pipeline$(tput sgr 0)"
        exit 0
    fi
    
    # TRYCYCLER CONTIGS CLUSTERING
    # running trycycler
    cd $BASEDIR/$CURRENT_DIR
    trycycler cluster --assemblies assemblies/*.fasta --reads longreads/${CURRENT_SUFFIX}.tmd.flt.fastq.gz --out_dir trycycler
    
    # quality checkpoint
    read -p "$(tput setaf 3)Please examine clusters and remove bad clusters. Continue? [y/n] $(tput sgr 0)" YN
    echo $YN
    if [ $YN != y ]
    then
        echo "$(tput setaf 3)Stopping pipeline$(tput sgr 0)"
        exit 0
    fi    

    # TRYCYCLER CLUSTER RECONCILING
    reconciling $i $t
}

i=36 # barcode name in '%02g'
t=19 # number of threads
j=200 # min read length

# pipeline $i $t $j

# for i in $(seq -f "%02g" 3 14)
# do
#       echo "$(tput setaf 3)Running pipeline with barcode$i$(tput sgr 0)"  
#       pipeline $i $t $j
# done

if [ -d "$BASEDIR/barcode$i/trycycler/cluster_001/1_contigs" ];
then
    echo "$(tput setaf 3)Cluster already exists, going to cluster reconciling with barcode0$i$(tput sgr 0)"
    reconciling $i $t
else
    echo "$(tput setaf 3)Running pipeline with barcode$i$(tput sgr 0)"
    pipeline $i $t $j
fi
