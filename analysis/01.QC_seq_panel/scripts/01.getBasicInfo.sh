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
      picard SamToFastq I=${bamOrg} FASTQ=${fastq%.gz}
   fi

   if [ ! -f $bamMem ]
   then
      bwa mem ${refGenome} ${fastq} | samtools view -buT ${refGenome} - | samtools sort - ${bamMem%.bam}
   fi

   if [ ! -f $bamAln ]
   then
      bwa aln ${refGenome} ${fastq} > ${bamAln%.bam}.sai
      bwa samse ${refGenome} ${bamAln%.bam}.sai ${fastq} | samtools view -buT ${refGenome} - | samtools sort - ${bamAln%.bam}
   fi
   echo $bamAln
done

exit 1

## Basic info
## ----------
if [ ! -f basic.info.txt ]
then
   for bam in ${bamPath}/*/*.bam
   do
      unmapped=`samtools view -f 4 $bam | wc -l`
      mapped=`samtools view -F 4 $bam | wc -l`
      echo -e "$bam\t$mapped\t$unmapped"
   done | awk '
      BEGIN{FS="\t";OFS="\t"}
      FNR==1{print "Bam","mapped_cnt","unmapped_cnt","mapped_perc","unmapped_perc"}
      {print $1,$2,$3,$2/($2+$3+1)*100,$3/($2+$3+1)*100}
      
   ' > basic.info.txt
fi

## index bam files
## ---------------
for bamFile in ${bamPath}/*.bam
do
   if [ ! -f ${bamFile}.bai ]
   then
      samtools index $bamFile
   fi
   name=`basename $bamFile`
   name=${name%.bam}

   ## Collect target Metrics
   if [ ! -d ${resultPath}/${name} ]
   then
      mkdir ${resultPath}/${name}
   fi

   ls -alh ${resultPath}/${name}/
   ls -alh ${resultPath}/${name}/${name}.10x.offtarget.bed

   if [ ! -f ${resultPath}/${name}/${name}.10x.offtarget.bed ]
   then

      ## Target Metric
      picard CollectTargetedPcrMetrics I=${bamFile} \
         O=${resultPath}/${name}/${name}.metrics.txt \
         R=${refGenome} \
         AMPLICON_INTERVALS=${targetFileList} \
         TARGET_INTERVALS=${targetFileList} \
         PER_TARGET_COVERAGE=${resultPath}/${name}/${name}.metrics.perTarget.txt

      ## Off-target metric
      bedtools genomecov -ibam ${bamFile} -g ${refGenome}.fai -bg | awk 'BEGIN{FS="\t";OFS="\t"} $4>9{print}' | bedtools merge -i - > ${resultPath}/${name}/${name}.10x.bed
      bedtools intersect -v -a ${resultPath}/${name}/${name}.10x.bed -b ${targetFileBed} > ${resultPath}/${name}/${name}.10x.offtarget.bed

      if [ `cat ${resultPath}/${name}/${name}.10x.offtarget.bed | wc -l` -gt 0 ]
      then
         picard BedToIntervalList I=${resultPath}/${name}/${name}.10x.offtarget.bed O=${resultPath}/${name}/${name}.10x.offtarget.list SD=${refGenomeDict}

         picard CollectTargetedPcrMetrics I=${bamFile} \
            O=${resultPath}/${name}/${name}.10x.offtarget.metrics.txt \
            R=${refGenome} \
            AMPLICON_INTERVALS=${resultPath}/${name}/${name}.10x.offtarget.list \
            TARGET_INTERVALS=${resultPath}/${name}/${name}.10x.offtarget.list \
            PER_TARGET_COVERAGE=${resultPath}/${name}/${name}.10x.offtarget.metrics.perTarget.txt
      else
         touch noOfftargets.txt
      fi
   fi
done




