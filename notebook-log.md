# Code for semester project
### Description of data
Source: https://elifesciences.org/articles/35962#content and https://www.ncbi.nlm.nih.gov/bioproject/PRJNA412631
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
Processed fastas retrieved from GitHub page:
https://github.com/lauringlab/Host_level_IAV_evolution

### Create concatenated fasta file for each segment
Created "concatenate.sh", which effectively grabs each segment from each sample and concatenates it to one file, per gene. Here's small sample:
```shell
sed -n '/>*NP/p' ./${f}.fasta > name
sed "s/NP/NP|${f}/g" name > name.temp
cat name.temp >> segmented_compiled-NP
sed -n '/>*NP/,/>/p' ./${f}.fasta > seq
sed '$ d' seq > seq.temp
sed '1d' seq.temp > seq
cat seq >> segmented_compiled-NP
echo -ne "[#--------] ${f}\r"
```

### Align HK_1 with ClustalW and Muscle
Created "run_clustalw.sh" as well as a muscle derivative following this format:
```shell
#!/bin/bash

## TIME START
SECONDS=0

#data_dir="/home/rieshunter/GitHub/myProject/data/aligned"
cwd=$(pwd)
echo $cwd

for file in segmented_compiled-*.fasta
do
f=${file%%.fasta}
echo $f
clustalw \
  -ALIGN \
  -INFILE=${cwd}/${f}.fasta \
  -OUTFILE=${cwd}/clustalw_${f}.fasta \
  -OUTPUT=FASTA
done

mkdir clustalw
mv clustalw* clustalw
mv *.dnd clustalw

## REPORT TIME             
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "----- Done -----"
echo ""
```

Benchmarking 10 flu genomes:
 - Muscle: 1 minute and 14 seconds
 - Clustalw: 3 minutes and 11 seconds
 - Tcoffee: [soon-to-be-clocked]
