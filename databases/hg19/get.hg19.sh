#!/bin/bash

## M Moisse - GenAdvice

if [ ! -f ucsc.hg19.1-22.X.Y.small.M.fa.gz ]
then
   for C in {1..22} X Y M;
   do
      if [[ "${C}" == "M" ]]
      then
         curl "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr${C}.fa.gz" | gunzip -c; 
      else 
         curl "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${C}.fa.gz" | gunzip -c; 
         a="d"
      fi
   done | bgzip > ucsc.hg19.1-22.X.Y.small.M.fa.gz
fi

if [ ! -f ucsc.hg19.1-22.X.Y.small.M.dict ]
then
   picard CreateSequenceDictionary R=ucsc.hg19.1-22.X.Y.small.M.fa.gz O=ucsc.hg19.1-22.X.Y.small.M.dict
fi

if [ ! -f ucsc.hg19.1-22.X.Y.small.M.fa.gz.fai ]
then
   samtools faidx ucsc.hg19.1-22.X.Y.small.M.fa.gz
fi

if [ ! -f ucsc.hg19.1-22.X.Y.small.M.fa.gz.pac ]
then
   bwa index ucsc.hg19.1-22.X.Y.small.M.fa.gz
fi

