#!/bin/bash

## GenAdvice - M Moisse
## Aug 2016

export GenAdviceHome=`pwd`
export PATH="${GenAdviceHome}/Ysoftware/anaconda/anaconda_3-4.1.1/anaconda3/bin/:$PATH"

cd analysis/01.QC_seq_panel
bash scripts/01.getBasicInfo.sh

cd $GenAdviceHome
