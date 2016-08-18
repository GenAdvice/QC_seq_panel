# QC_seq_panel

## Installation
The installation requires 9Gib of disk space and 8Gib of memory and takes 1 to 2 hours to finish.
```bash
git clone https://github.com/GenAdvice/QC_seq_panel.git
cd QC_seq_panel
bash init.env.sh
```

## Running the tool
Add the bam files you want to analyse in QC_seq_panel/analysis/01.QC_seq_panel/data/bam/org/
Also add the target file in QC_seq_panel/analysis/01.QC_seq_panel/data/ and call it targetFile.bed

Running the analysis takes about 5min per sample
```bash
## in the QC_seq_panel directory
bash run.analysis.sh
`
