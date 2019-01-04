#!/bin/bash
# Deconvolution project data
# Process fastq data by running cutadapt, STAR, featureCounts 

# software indicator
Run_cutadapt=0
Run_STAR=0
Run_featureCounts=1


#==============================================================================================
# parameter for cutadapt
#==============================================================================================
FASTQDIR=/ddn/gs1/group/li3/DataFromLLI/DeconvData/NS50593_94_NS50596_97 
FASTQTRIMMED=/ddn/gs1/group/li3/DataFromLLI/DeconvData/NS50593_94_NS50596_97_trimmed  

#==============================================================================================
# parameters for STAR
#==============================================================================================
OUTFILTERMISMATCHNOVERLMAX=0.04
RUNTHREADN=14
RUNMODE=alignReads
GENOMEDIR=/ddn/gs1/group/li3/DataFromLLI/DeconvData/GTF/
OUTSAMTYPE_1=BAM
OUTSAMTYPE_2=Unsorted
OUTSAMSTRANDFIELD=intronMotif
OUTMULTIMAPPERORDER=Random
OUTSAMATTRIHSTART=0
OUTFILTERTYPE=BySJout
ALIGNSJOVERHANMIN=8
QUANTMODE=GeneCounts
READFILESIN=$FASTQTRIMMED
STAR_OUTPUT=/ddn/gs1/group/li3/DataFromLLI/DeconvData/STAR_output_trimmed


#==============================================================================================
# parameters for featureCounts
#==============================================================================================
ANNOTATION=/ddn/gs1/group/li3/DataFromLLI/DeconvData/GTF/Homo_sapiens.GRCh38.90.gtf
ANNOTATIONHG=/ddn/gs1/group/li3/DataFromLLI/DeconvData/GTF/hg19_RefSeq_exon.txt  
NTHREADS=14
FEATURETYPE=exon
#ATTRIBUTE=gene_id
ATTRIBUTE=gene_name
BAMDIR=$STAR_OUTPUT/STAR_BAM
FCDIR=/ddn/gs1/group/li3/DataFromLLI/DeconvData/FeatureCounts_output_trimmed
FCGENENAMEDIR=/ddn/gs1/group/li3/DataFromLLI/DeconvData/FeatureCounts_output_trimmed_gene_name

#==============================================================================================
# call cutadapt
#==============================================================================================
if ((Run_cutadapt==1)) 
then
  printf "Running cutadapt...\n"
  for ((filename=1;filename<=40;filename++))
  {
    cutadapt -o $FASTQTRIMMED/IS-$filename-trimmed.fastq -q 20 -m 18 -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC $FASTQDIR/IS-$filename.fastq
  }
fi

#==============================================================================================
# call STAR
#==============================================================================================
if ((Run_STAR==1)) 
then
  printf "Running STAR...\n"
  for ((filename=1;filename<=40;filename++))
  {
    STAR --genomeDir $GENOMEDIR --readFilesIn $READFILESIN/IS-$filename-trimmed.fastq --outFilterMismatchNoverLmax $OUTFILTERMISMATCHNOVERLMAX --runThreadN $RUNTHREADN --outSAMtype $OUTSAMTYPE_1 $OUTSAMTYPE_2 --outFileNamePrefix $STAR_OUTPUT/IS-$filename-trimmed.fastq
  }
  mv $STAR_OUTPUT/*.bam $STAR_OUTPUT/STAR_BAM
fi

#==============================================================================================
# call featureCounts 
#==============================================================================================
if ((Run_featureCounts==1))
then
  printf "Running featureCounts...\n"
  for ((filename=1;filename<=40;filename++))
  {
    if (($ATTRIBUTE=="gene_id"))
    then  
    printf "Using gene_id as attribute..."
    featureCounts -a $ANNOTATION -T $NTHREADS -t $FEATURETYPE -g $ATTRIBUTE -o $FCDIR/IS-$filename-trimmed.fastq.counts.txt $BAMDIR/IS-$filename-trimmed.fastqAligned.out.bam
    fi

    if (($ATTRIBUTE=="gene_name"))
    then
    printf "Using gene_name as attribute..."
    featureCounts -a $ANNOTATION -T $NTHREADS -t $FEATURETYPE -g $ATTRIBUTE -o $FCGENENAMEDIR/IS-$filename-trimmed-gene-name.fastq.counts.txt $BAMDIR/IS-$filename-trimmed.fastqAligned.out.bam
    fi
  }
fi 


