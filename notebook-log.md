# Code for semester project
### Description of data
Source: https://elifesciences.org/articles/35962#content
Who: HIVE Cohort for respiratory virus surveillance
What: Flu sequences from 200 people in a test-negative study
Where: JT McCrone et al. at UMich
When: 2010-2015
Why: We know flu evolution is evolutionarily stochastic within-host but deterministic globally. Maybe the answer lies between-host?

Aims:
Evolution w/in hosts through serial time points
Evolution b/w hosts through household pairing

### Stage of data
Raw SRA reads in need of QC and alignment

### Download data
Downloaded SRA toolkit
Created SraAccList.csv from NCBI source
Created short version with:
```shell
cat SraAccList.csv | head -6 | tail -5 > short_SraAccList.csv
```
Download the first five SRA files:
```shell
prefetch --option-file short_SraAccList.csv
```
Convert to fastq, organize, compress:
```shell
fasterq-dump --split-files ./SRR*/*.sra
rm -r SRR*/
mkdir ./raw_reads
mv SRR*.fastq ./raw_reads
cd ./raw_reads
gzip *.fastq #this part takes a while
```

### Trim reads
```shell
TrimmomaticPE


```








