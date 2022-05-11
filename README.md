

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

Output for mummer is <prefix>.contigs.fasta (The genome)
  
## FASTA QC
```
module load bioinfo-tools samtools/0.1.19 bwa
module load FastQC

# Quality Before
fastqc /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/genomics_data/Illumina/* -o /home/samuel50/Genom_Analysis/Genomic_Data/Illumina/FASTAQC/Quality_before

# Quality After
fastqc /home/samuel50/Genom_Analysis/Genomic_Data/Illumina/Trimmomatic/* -o /home/samuel50/Genom_Analysis/Genomic_Data/Illumina/FASTAQC/Quality_after

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
## BLAST

### Command for blasting sequence with latest (2.9.0+) nucleotide database.
 ```
 #!/bin/bash
#SBATCH -A uppmax2022-2-5
#SBATCH -M snowy
#SBATCH -p core -n 2
#SBATCH -t 08:00:00
#SBATCH -J BLAST-04-09


#load modules
module load bioinfo-tools samtools/0.1.19 bwa
module load blast/2.9.0+

#Command
blastn -db nt -query Genome_assembly_Pacbio-30-03-16.38.contigs.fasta -out out_BLAST.txt
 ```

From the output we can see that "Enterococcus faecium strain E745, complete genome" has the best alignment to our sequence from the PAC-BIO output. Thus this is Enterococcus faecium strain E745, GenBank: CP014529.1 is your reference moving forward with analysis.
 
 
 
## Mummer and Mummerplot
Before running mummerplot the user needs to align a reference with the contigs from the Genome assembly (canu) with mummer.


```
module load bioinfo-tools samtools/0.1.19 bwa
module load MUMmer
mummer -mum -b -c /domus/h1/samuel50/Genom_Analysis/Genomic_Data/Reference/sequence.fasta /domus/h1/samuel50/Genom_Analysis/Genomic_Data/PACBIO/Complete_Assembly/Complete_Assembly_PACKBIO_30-03-16.38/Genome_assembly_Pacbio-30-03-16.38.contigs.fasta -save mum.txt > ref_qry.mums
 
 mummerplot --postscript --prefix=ref_qry ref_qry.mums

gnuplot ref_qry.gp

```
Then view or print the postscript plot ref_qry.ps in whatever manner you wish.

 ![ref_qry](https://user-images.githubusercontent.com/370074/163407531-283b7eec-4d74-467f-884a-f4eeba89a62d.png)

 
## Prokka
 
 ```
 #!/bin/bash
#SBATCH -A uppmax2022-2-5
#SBATCH -M snowy
#SBATCH -p core -n 2
#SBATCH -t 08:00:00
#SBATCH -J Prokka-2
#SBACTH --reservation=uppmax2022-2-5_8

#load modules
module load bioinfo-tools samtools/0.1.19 bwa prokka

prokka \
--outdir /home/samuel50/Genom_Analysis/Prokka/result \
--addgenes --cpus 2 --kingdom Bacteria --species Enterococcus_faecium --strain E745 --prefix mygenome \
/domus/h1/samuel50/Genom_Analysis/Genomic_Data/PACBIO/Complete_Assembly/Complete_Assembly_PACKBIO_30-03-16.38/Genome_assembly_Pacbio-30-03-16.38.contigs.fasta

```
 
## BAM
 ```
 
#!/bin/bash
#SBATCH -A uppmax2022-2-5
#SBATCH -M snowy
#SBATCH -p core -n 2
#SBATCH -t 08:00:00
#SBATCH -J BWA


#load modules
module load bioinfo-tools samtools/0.1.19 bwa


#Command BH and Serum

bwa mem /domus/h1/samuel50/Genom_Analysis/Genomic_Data/PACBIO/Complete_Assembly/Complete_Assembly_PACKBIO_30-03-16.38/Genome_assembly_Pacbio-30-03-16.38.contigs.fasta\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797972_pass_1.fastq.gz\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797972_pass_2.fastq.gz | samtools view -bS - > /domus/h1/samuel50/Genom_Analysis/BWA/BH-ERR1797972-file.$

bwa mem /domus/h1/samuel50/Genom_Analysis/Genomic_Data/PACBIO/Complete_Assembly/Complete_Assembly_PACKBIO_30-03-16.38/Genome_assembly_Pacbio-30-03-16.38.contigs.fasta\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797973_pass_1.fastq.gz\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797973_pass_2.fastq.gz | samtools view -bS - > /domus/h1/samuel50/Genom_Analysis/BWA/BH-ERR1797973-file.$

bwa mem /domus/h1/samuel50/Genom_Analysis/Genomic_Data/PACBIO/Complete_Assembly/Complete_Assembly_PACKBIO_30-03-16.38/Genome_assembly_Pacbio-30-03-16.38.contigs.fasta\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797974_pass_1.fastq.gz\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797974_pass_2.fastq.gz | samtools view -bS - > /domus/h1/samuel50/Genom_Analysis/BWA/BH-ERR1797974-file.$



bwa mem /domus/h1/samuel50/Genom_Analysis/Genomic_Data/PACBIO/Complete_Assembly/Complete_Assembly_PACKBIO_30-03-16.38/Genome_assembly_Pacbio-30-03-16.38.contigs.fasta\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797969_pass_1.fastq.gz\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797969_pass_2.fastq.gz | samtools view -bS - > /domus/h1/samuel50/Genom_Analysis/BWA/Serum-ERR1797969$

bwa mem /domus/h1/samuel50/Genom_Analysis/Genomic_Data/PACBIO/Complete_Assembly/Complete_Assembly_PACKBIO_30-03-16.38/Genome_assembly_Pacbio-30-03-16.38.contigs.fasta\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797970_pass_1.fastq.gz\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797970_pass_2.fastq.gz | samtools view -bS - > /domus/h1/samuel50/Genom_Analysis/BWA/Serum-ERR1797970$

bwa mem /domus/h1/samuel50/Genom_Analysis/Genomic_Data/PACBIO/Complete_Assembly/Complete_Assembly_PACKBIO_30-03-16.38/Genome_assembly_Pacbio-30-03-16.38.contigs.fasta\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797971_pass_1.fastq.gz\
 /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797971_pass_2.fastq.gz | samtools view -bS - > /domus/h1/samuel50/Genom_Analysis/BWA/Serum-ERR1797971$



   
hexdump -C BH-ERR1797972-file.bam 
 
 ```
 
 
 ## HTseq
 
 ```
#!/bin/bash
#SBATCH -A uppmax2022-2-5
#SBATCH -M snowy
#SBATCH -p core -n 2
#SBATCH -t 06:00:00
#SBATCH -J HT-seq_Serum


#load modules
module load bioinfo-tools samtools/0.1.19 bwa htseq

#location to BAM files 
cd /home/samuel50/Genom_Analysis/BWA

htseq-count -f bam -r pos -i ID -s no -t CDS Serum-ERR1797969-file.bam Serum-ERR1797970-file.bam Serum-ERR1797971-file.bam /home/samuel50/Genom_Analysis/Prokka/result/mygenome.gff> counts-Serum.txt
 
 
 htseq-count -f bam -r pos -i ID -s no -t CDS BH-ERR1797972-file.bam BH-ERR1797973-file.bam BH-ERR1797974-file.bam /home/samuel50/Genom_Analysis/Prokka/result/mygenome.gff> counts-BH.txt
 
 ```
 
## Cut my files
 
 ```
 #filed 1 and 4 for example
 cut -d$'\t' -f 1,4 counts-Serum.txt >counts-Serum-3.txt
 
 ```
 
## DEseq
```
library("DESeq2")
directory <- "/home/samuel50/Genom_Analysis/htseq"

sampleFiles <- c('counts-BH-1.txt','counts-BH-2.txt','counts-BH-3.txt','counts-Serum-1.txt','counts-Serum-2.txt','counts-Serum-3.txt')
sampleCondition <- c('BH','BH','BH','Serum','Serum','Serum')
sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = sampleCondition)

sampleTable$condition <- factor(sampleTable$condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq

dds <- DESeq(ddsHTSeq)

res <- results(dds)

summary(res)

plotMA(res)

library(ReportingTools)
report <- HTMLReport(shortName = 'Differential expression analysis BH vs Serum', title = 'Differential expression analysis BH vs Serum', reportDirectory = '.')
publish(dds, report, pvalueCutoff=0.05, make.plots = TRUE, factor = sampleTable$condition, reportDir = ".")
finish(report) 
```
![Untitled](https://user-images.githubusercontent.com/370074/167854184-90cb323f-bbe9-44b7-91e3-0f70a69a37af.png)

