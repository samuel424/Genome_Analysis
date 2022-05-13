
README

  This github contains the bioinfromatic scrips and tools inorder to reprduce the pipline of the article
  "RNA-seq and Tn-seq reveal fitness determinants of vancomycin-resistant Enterococcus faecium during growth in human serum". 

PIPLINE

  This is the order of the Scripts for the Analysis, and to how use them:

  1. Genome assemmbly with Canu (Pac-Bio data)
  
  code: sbatch sbatch_Assembly.

  2. Reads processing with Trimmomatic (Illumina data)
  
  code: sbatch sbatch_Trimmo.

  3. Check quality with FASTA QC for the Trimmomatic data (before and after)
  
  code: sbatch_FASTAQC

  4. Use BLAST to find genome
  
  code: sbatch sbatch_BLAST

  6. Mummer and Mummerplot for assembly evaluation
  
  Before running mummerplot the user needs to align a reference with the contigs from the Genome assembly (canu) with mummer.
  
  code: sbatch sbatch_Mummer

  8. Artemis Comparison Tool (ACT) is another visualization tool 
  
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

  cut -d$'\t' -f {fileds} counts-{condition}.txt > {output}txt

  14 Use DEseq to compute a differential expression analysis between rich medium an heat-inectivated serum codntion
  
  code: sbatch sbatch_DEseq

CONTACT

  If you have problems or questions please contact me at samuel.50@live.se
  
GIT

  to download the very latest cource off the GIT server do this: 

  git clone https://github.com/samuel424/Genom_Analysis.git

```
  13. The file created in from the HTseq sbatch job need to be cut wiiht the folloing code
