#!/bin/bash
# Deconvolution project data
# Run cat or zcat to combine two fastq files that come from the same library


FASTQDIR_1=/ddn/gs1/project/nextgen/post/li3/NS50593-Li-Shats
FASTQDIR_2=/ddn/gs1/project/nextgen/post/li3/NS50594-Li-Shats
FASTQDIR_3=/ddn/gs1/project/nextgen/post/li3/NS50596-Li-Shats
FASTQDIR_4=/ddn/gs1/project/nextgen/post/li3/NS50597-Li-Shats

FASTQOUT=/ddn/gs1/group/li3/DataFromLLI/DeconvData/NS50593_94_NS50596_97 

#zcat $FASTQFILE_11 $FASTFILE_21 | gzip -c > IS-1-zcat.fastq.gz

for ((filename=1;filename<=40;filename++))
{
if (($filename <=20)) 
then
  FASTQFILE_1=$FASTQDIR_1/NS50593_171103_NS500489_AHNLM2BGX3.IS-$filename.single.sanger.fastq.gz
  FASTQFILE_2=$FASTQDIR_2/NS50594_171106_NS500489_AHNLHTBGX3.IS-$filename.single.sanger.fastq.gz
  #printf "$FASTQFILE_1\n$FASTQFILE_2\n\n"
  cat $FASTQFILE_1 $FASTFILE_2 > $FASTQOUT/IS-$filename.fastq.gz
else
  FASTQFILE_3=$FASTQDIR_3/NS50596_171108_NS500489_AHNMYFBGX3.IS-$filename.single.sanger.fastq.gz 
  FASTQFILE_4=$FASTQDIR_4/NS50597_171109_NS500489_AHMWTHBGX3.IS-$filename.single.sanger.fastq.gz
  #printf "$FASTQFILE_3\n$FASTQFILE_4\n\n"
  cat $FASTQFILE_3 $FASTFILE_4 > $FASTQOUT/IS-$filename.fastq.gz
fi
}

exit 0
