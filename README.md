# Genom_Analysis

### Clone Github
Inorder to clone github use  git clone https://github.com/samuel424/Genom_Analysis.git

This is example for git add,commit and push for test1_30_03_2022 in the map Script
```
git add test1_30-03-2022

git commit -m "added testfile"

git push
```
If the user wants to pull data from Github use the code ```git pull```

## Assemmbly Pac-BIO
```
#!/bin/bash
#SBATCH -A uppmax2022-2-5
#SBATCH -M snowy
#SBATCH -p core
#SBTACH -n 2
#SBATCH -t 07:00:00
#SBATCH -J Genome_Assembly_PacBIO_2022-03-30-16:38
#SBATCH -o Genome_Assembly_PacBIO-2022-03-30.log 
#SBATCH -e Genome_Assembly_PacBIO-2022-03-30.err 

#load modules
module load bioinfo-tools samtools/0.1.19 bwa
module load canu/2.0

#Command

canu \
 -p Genome_assembly_Pacbio-30-03-16.38 -d /home/samuel50/Genom_Analysis/Genomic_Data  \
 genomeSize=2.4m \
 -pacbio /home/samuel50/Genom_Analysis/Genomic_Data/PACBIO_data/data/*
```
## FASTA QC
```
module load bioinfo-tools samtools/0.1.19 bwa
module load canu/2.0
fastaqc /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/genomics_data/Illumina/* -o /home/samuel50/Genom_Analysis/Genomic_Data/Illumina/FASTAQC
```
## Trimmomatic
```
#!/bin/bash
#SBATCH -A uppmax2022-2-5
#SBATCH -M snowy
#SBATCH -p core
#SBTACH -n 2
#SBATCH -t 02:00:00
#SBATCH -J Genome_Trimmo-2022-03-30-15:15
#SBATCH -o Genome_Trimmo-2022-03-30.log
#SBATCH -e Genome_Trimmo-2022-03-30.err

#load modules
module load bioinfo-tools samtools/0.1.19 bwa
module load trimmomatic

#Command
 
trimmomatic PE -phred33 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz
/proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz paired.fq.gz_p1 fq.gz_s1 fq.gz_p2 
fq.gz_s2 ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

```
![Screenshot 2022-05-11 at 14 02 20](https://user-images.githubusercontent.com/370074/167852586-04161b12-15c3-4a1b-bf7f-5764caac2348.png)

