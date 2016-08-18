#!/bin/bash

## GenAdvice - M Moisse
## Aug 2016

export GenAdviceHome=`pwd`

## Get the right software
bash software/getSoftware.sh

export PATH="${GenAdviceHome}/software/anaconda/anaconda_3-4.1.1/anaconda3/bin/:$PATH"

## Get the right ref genome
bash databases/hg19/get.hg19.sh
