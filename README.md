```
README

This github contains the bioinfromatic scrips and tools inorder to reprduce the pipline of the article "RNA-seq and Tn-seq reveal fitness determinants of vancomycin-resistant Enterococcus faecium during growth in human serum". 

This is the order of the Scripts for the Analysis, and to how use them:

1.Genome assemmbly with Canu (Pac-Bio data)
Start with the DNA assenmbly with the sbatch job named sbatch_Assembly
code: sbatch sbatch_Assembly.

2.Reads processing with Trimmomatic (Illumina data)
code: sbatch sbatch_Trimmo.

3. Check quality with FASTA QC for the Trimmomatic data (before and after) with the following command :


module load bioinfo-tools samtools/0.1.19 bwa
module load FastQC

# Quality Before
fastqc /proj/uppmax2022-2-5/Genome_Analysis/1_Zhang_2017/genomics_data/Illumina/* -o /home/samuel50/Genom_Analysis/Genomic_Data/Illumina/FASTAQC/Quality_before

# Quality After
fastqc /home/samuel50/Genom_Analysis/Genomic_Data/Illumina/Trimmomatic/* -o /home/samuel50/Genom_Analysis/Genomic_Data/Illumina/FASTAQC/Quality_after

4. Use BLAST to find genome
code: sbatch sbatch_BLAST

6.  Mummer and Mummerplot for assembly evaluation
Before running mummerplot the user needs to align a reference with the contigs from the Genome assembly (canu) with mummer.
Use the folloing code in comandline:

module load bioinfo-tools samtools/0.1.19 bwa
module load MUMmer
mummer -mum -b -c /domus/h1/samuel50/Genom_Analysis/Genomic_Data/Reference/sequence.fasta /domus/h1/samuel50/Genom_Analysis/Genomic_Data/PACBIO/Complete_Assembly/Complete_Assembly_PACKBIO_30-03-16.38/Genome_assembly_Pacbio-30-03-16.38.contigs.fasta -save mum.txt > ref_qry.mums
 
 mummerplot --postscript --prefix=ref_qry ref_qry.mums

gnuplot ref_qry.gp



8. Artemis Comparison Tool (ACT) is another visualization tool especially designed for displaying pairwise comparisons between two or more DNA sequences
If you want to run Artemis/ACT from Uppmax, you need to keep this into mind:
-Use ssh -AX to log into Uppmax
-Use “interactive” instead of “salloc” to allocate cores
-Don’t use the Windows command line, try instead an interface like MobaXTerm or PuTTy

9. Prokka for structural and functional annotation
code: sbatch sbatch_Prokka 

11. BWA to prefomre the RNA-Seq reads aligmment against the assmebled geome
code : sbatch sbatch_BWA

12. Use HTseq to count the amount of reads against the whole gneome
code: sbatch sbatch_HTseq

13. The file created in from the HTseq sbatch job need to be cut wiiht the folloing code
Filed 1 and 4 for example below:
cut -d$'\t' -f 1,4 counts-Serum.txt >counts-Serum-3.txt

14 Use DEseq to compute a differential expression analysis between rich medium an heat-inectivated serum codntion
code: sbatch sbatch_DEseq


```
