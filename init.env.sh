#!/bin/bash

## GenAdvice - M Moisse
## Aug 2016

export GenAdviceHome=`pwd`

## Get the right software
cd software
bash getSoftware.sh
cd $GenAdviceHome

export PATH="${GenAdviceHome}/software/anaconda/anaconda_3-4.1.1/anaconda3/bin/:$PATH"

## Get the right ref genome
cd databases/hg19/
bash get.hg19.sh
cd $GenAdviceHome
