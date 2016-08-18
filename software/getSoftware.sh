#!/bin/bash

softHome=`pwd`

## install aconda
tool=anaconda
vernbr=3-4.1.1
version=${tool}_${vernbr}

cd ${softHome}

if [ ! -d ${tool} ]
then
   mkdir ${tool}
fi

cd ${tool}

if [ ! -d ${version} ]
then
   mkdir ${version}
   cd ${version}
   #wget http://repo.continuum.io/archive/Anaconda3-4.1.1-Linux-x86_64.sh

   bash Anaconda3-4.1.1-Linux-x86_64.sh -b -p ${GenAdviceHome}/software/${tool}/${version}/anaconda3

   # . ~/.bashrc
   export PATH="${GenAdviceHome}/software/${tool}/${version}/anaconda3/bin:$PATH"
   echo $PATH

   conda config --add channels r
   conda config --add channels BioBuilds
   conda config --add channels bioconda

fi


## picard tools
conda install -c bioconda picard

## samtools
conda install -c bioconda samtools

## igv
conda install igv

## bwa
conda install bwa

## bwa
conda install bedtools

