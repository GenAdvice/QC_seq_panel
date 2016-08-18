#!/bin/bash

## M Moisse - GenAdvice
## Aug 2016

## =========
## Variables
## =========
refGenome=${GenAdviceHome}/databases/hg19/ucsc.hg19.1-22.X.Y.small.M.fa.gz
refGenomeDict=${GenAdviceHome}/databases/hg19/ucsc.hg19.1-22.X.Y.small.M.dict
targetFileBed=data/targetFile.bed
targetFileList=data/targetFile.list

if [ ! -f ${targetFileBed} ]
then
   echo "Please add the target file! `pwd`/${targetFileBed}"
   exit 1
fi

if [ $# -ge 2 ]
then
   bamPath=$1
   resultPath=$2
else
   bamPath=data/bam
   resultPath=results
fi

## Convert target bed-file to list-file
picard BedToIntervalList I=${targetFileBed} O=${targetFileList} SD=${refGenomeDict}

## ====
## Code 
## ====

if [ ! -d ${bamPath}/mem/ ]
then
   mkdir ${bamPath}/mem/
fi

if [ ! -d ${bamPath}/aln/ ]
then
   mkdir ${bamPath}/aln/
fi

if [ ! -d ${bamPath}/fastq/ ]
then
   mkdir ${bamPath}/fastq/
fi

if [ ! -d ${resultPath} ]
then
   mkdir ${resultPath}
fi

## Remap bam files
## ---------------
for bamOrg in ${bamPath}/org/*.bam
do
   bamMem=${bamOrg/org/mem}
   bamAln=${bamOrg/org/aln}
   fastq=${bamOrg/org/fastq}
   fastq=${fastq%.bam}.fastq.gz

   if [ ! -f ${fastq} ]
   then
      echo "creating fastq from bam"
      picard SamToFastq I=${bamOrg} FASTQ=${fastq%.gz}
      gzip ${fastq%.gz}
   fi

   if [ ! -f $bamMem ]
   then
      echo "mapping with bwa mem"
      bwa mem ${refGenome} ${fastq} | samtools view -buT ${refGenome} - | samtools sort - > ${bamMem}
   fi

   if [ ! -f $bamAln ]
   then
      echo "mapping with bwa aln"
      bwa aln ${refGenome} ${fastq} > ${bamAln%.bam}.sai
      bwa samse ${refGenome} ${bamAln%.bam}.sai ${fastq} | samtools view -buT ${refGenome} - | samtools sort - > ${bamAln}
   fi
   echo $bamAln
done


## Basic info
## ----------

## TODO only calc that of new bam files
# if [ ! -f ${resultPath}/basic.info.txt ]
# then
   echo "Getting basic info"
   for bam in ${bamPath}/*/*.bam
   do
      # x=`grep aaaaa targetFile.list` 
      # if [ -z "$x" ]; then echo "no"; fi

      unmapped=`samtools view -f 4 $bam | wc -l`
      mapped=`samtools view -F 4 $bam | wc -l`
      echo -e "$bam\t$mapped\t$unmapped"
   done | awk '
      BEGIN{FS="\t";OFS="\t"}
      FNR==1{print "Bam","mapped_cnt","unmapped_cnt","mapped_perc","unmapped_perc"}
      {print $1,$2,$3,$2/($2+$3+1)*100,$3/($2+$3+1)*100}
      
   ' > ${resultPath}/basic.info.txt
# fi

## index bam files
## ---------------
for bamFile in ${bamPath}/*/*.bam
do
   if [ ! -f ${bamFile}.bai ]
   then
      samtools index $bamFile
   fi
   name=`basename $bamFile`
   name=${name%.bam}

   mapping=`dirname $bamFile`
   mapping=`basename $mapping`


   ## Collect target Metrics
   ## ----------------------
   workPath=${resultPath}/${mapping}
   if [ ! -d ${workPath} ]
   then
      mkdir ${workPath}
   fi
   
   workPath=${workPath}/${name}
   if [ ! -d ${workPath} ]
   then
      mkdir ${workPath}
   fi

   if [ ! -f ${workPath}/${name}.10x.offtarget.bed ]
   then

      ## Target Metric
      picard CollectTargetedPcrMetrics I=${bamFile} \
         O=${workPath}/${name}.metrics.txt \
         R=${refGenome} \
         AMPLICON_INTERVALS=${targetFileList} \
         TARGET_INTERVALS=${targetFileList} \
         PER_TARGET_COVERAGE=${workPath}/${name}.metrics.perTarget.txt

      ## Off-target metric
      bedtools genomecov -ibam ${bamFile} -g ${refGenome}.fai -bg | awk 'BEGIN{FS="\t";OFS="\t"} $4>9{print}' | bedtools merge -i - > ${workPath}/${name}.10x.bed
      bedtools intersect -v -a ${workPath}/${name}.10x.bed -b ${targetFileBed} > ${workPath}/${name}.10x.offtarget.bed

      if [ `cat ${workPath}/${name}.10x.offtarget.bed | wc -l` -gt 0 ]
      then
         picard BedToIntervalList I=${workPath}/${name}.10x.offtarget.bed O=${workPath}/${name}.10x.offtarget.list SD=${refGenomeDict}

         picard CollectTargetedPcrMetrics I=${bamFile} \
            O=${workPath}/${name}.10x.offtarget.metrics.txt \
            R=${refGenome} \
            AMPLICON_INTERVALS=${workPath}/${name}.10x.offtarget.list \
            TARGET_INTERVALS=${workPath}/${name}.10x.offtarget.list \
            PER_TARGET_COVERAGE=${workPath}/${name}.10x.offtarget.metrics.perTarget.txt
      else
         touch ${workPath}/noOfftargets.txt
      fi
   fi
done




