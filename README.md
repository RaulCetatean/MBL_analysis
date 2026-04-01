# MBL_analysis
Code I used for the bioinformatic analysis of the MBL isolates

## 1. SNAKEFILE

The Snakefile includes the commands run for each isolate starting from the raw reads.
The analysis goes as follows:
1. Quality control with FastQC
2. Adapter removal and quality trimming with trim_galore. As parameters used, minimum length of 50 and minimum Phred score of 20
3. FastQC on the trimmed reads
4. Kraken and Bracken for contamination check
5. De novo assembly with SPAdes
6. QUAST evaluation metrics
7. MLST
8. Kleborate for KL locus and O antigen, as well as AMR and virulence gene detection
9. Annotation with Bakta (in my case it has its own environment with only the Bakta program)
10. AMR gene detection with AMRFinderPlus
11. PlasmidFinder for in silico plasmid detection

## 2. filter_pathogenwatch_data.py

This code was used to extract only the genomes having co-occurrence of NDM-1 or NDM-5 with OXA-48

## 3. extract_iron_gene_sequences.py

extract_iron_gene_sequences.py extracts the amino acid sequences of the iron uptake/transport genes to obtain a multi-fasta file with the sequences of all isolates. The output will then be aligned with MAFFT to identify potential mutations

## 4. Create_heatmap_poster.R

create_heatmap_poster.R creates a heatmap on the data obtained by AMRFinderPlus. The data should be in .csv format and binary (1 or 0 for presence or absence).

# CITATIONS
[FASTQC]: https://github.com/s-andrews/FastQC
[TRIM_GALORE]: https://github.com/FelixKrueger/TrimGalore
[KRAKEN] Wood, D.E.; Lu, J.; Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol. 2019, 20, 257.
[BRACKEN]  Lu, J.; Breitwieser, F.P.; Thielen, P.; Salzberg, S.L. Bracken: estimating species abundance in metagenomics data. PeerJ Comput Sci. 2017, 3, e104.
[SPADES] Prjibelski, A.; Antipov, D.; Meleshko, D.; Lapidus, A.; Korobeynikov, A. Using SPAdes De Novo Assembler. Curr Protoc Bioinformatics. 2020, 70, e102. 
[QUAST] Mikheenko, A.; Saveliev, V.; Hirsch, P.; Gurevich, A. WebQUAST: online evaluation of genome assemblies. Nucleic Acids Res. 2023,51, W601-W606. 
[BAKTA] Schwengers, O.; Jelonek, L.; Dieckmann, M.A.; Beyvers, S.; Blom, J.; Goesmann, A. Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microb Genom. 2021, 7, 000685.
[MLST] Jolley, K.A.; Maiden, M.C. BIGSdb: Scalable analysis of bacterial genome variation at the population level. BMC Bioinform. 2010, 11, 595. 
[KLEBORATE1] Wyres, K.L.; Wick, R.R.; Gorrie, C.; Jenney, A.; Follador, R.; Thomson, N.R.; Holt, K.E. Identification of Klebsiella capsule synthesis loci from whole genome data. Microb Genom. 2016, 2, e000102. 
[KLEBORATE2] Lam, M.M.C.; Wick, R.R.; Watts, S.C.; Cerdeira, L.T.; Wyres, K.L.; Holt, K.E. A genomic surveillance framework and genotyping tool for Klebsiella pneumoniae and its related species complex. Nat Commun. 2021, 12, 4188.
[AMRFINDER] Feldgarden, M.; Brover, V.; Gonzalez-Escalona, N.; Frye, J.G.; Haendiges, J.; Haft, D.H.; Hoffmann, M.; Pettengill, J.B.; Prasad, A.B.; Tillman, G.E.; et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci. Rep. 2021, 11, 12728.
[PLASMIDFINDER] Carattoli, A.; Hasman, H. PlasmidFinder and In Silico pMLST: Identification and Typing of Plasmid Replicons in Whole-Genome Sequencing (WGS). Methods Mol. Biol. 2020, 2075, 285–294.
[SNAKEMAKE] Mölder, F.; Jablonski, K.P.; Letcher, B.; Hall, M.B.; Tomkins-Tinch, C.H.; Sochat, V.; Forster, J.; Lee, S.; Twardziok, S.O.; Kanitz, A.; et al. Sustainable data analysis with Snakemake. F1000Research 2021, 10, 33.
[PARSNP] Kille, B; Nute, M.G.; Huang, V.; Kim, E.; Phillippy, A.M.; Treangen, T.J. Parsnp 2.0: scalable core-genome alignment for massive microbial datasets. Bioinformatics. 2024, ;40, 1-5.  
[GUBBINS] Croucher, N.J.; Page, A.J.; Connor, T.R.; Delaney, A.J.; Keane, J.A.; Bentley, S.D.; Parkhill, J.; Harris, S.R. Rapid phylogenetic analysis of large samples of recombinant bacterial whole genome sequences using Gubbins. Nucleic Acids Res. 2015, 43, e15. 
[MAFFT] Katoh, K.; Standley, D.M. MAFFT multiple sequence alignment software version 7: Improvements in performance and usability. Mol. Biol. Evol. 2013, 30, 772–780.
