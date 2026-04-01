configfile: 'config/config.yaml'

rule all:
    input:
        #FASTQC
        expand('out/fastqc/raw/{sample}/{sample}_R1_001_fastqc.html',sample=config['samples']),
        expand('out/fastqc/raw/{sample}/{sample}_R2_001_fastqc.html',sample=config['samples']),
        expand('out/fastqc/raw/{sample}/{sample}_R1_001_fastqc.zip',sample=config['samples']),
        expand('out/fastqc/raw/{sample}/{sample}_R2_001_fastqc.zip',sample=config['samples']),

        # TRIM_GALORE
        expand("out/trim/{sample}_val_1.fq.gz",sample=config["samples"]),
        expand("out/trim/{sample}_val_2.fq.gz",sample=config["samples"]),

        #FASTQC ON TRIMMED READS
        expand("out/fastqc/trimmed/{sample}/{sample}_val_1_fastqc.html",sample=config["samples"]),
        expand("out/fastqc/trimmed/{sample}/{sample}_val_2_fastqc.html",sample=config["samples"]),
        expand("out/fastqc/trimmed/{sample}/{sample}_val_1_fastqc.zip",sample=config["samples"]),
        expand("out/fastqc/trimmed/{sample}/{sample}_val_2_fastqc.zip",sample=config["samples"]),

        #GUNZIP READS
        expand("data_unzipped/{sample}_L001_R1_001.fastq",sample=config["samples"]),
        expand("data_unzipped/{sample}_L001_R2_001.fastq",sample=config["samples"]),
        expand("data_unzipped/{sample}_val_1.fq",sample=config["samples"]),
        expand("data_unzipped/{sample}_val_2.fq",sample=config["samples"]),

        #CONTAMINATION CHECK: KRAKEN AND BRACKEN
        expand("out/kraken2/{sample}/{sample}.kreport",sample=config["samples"]),
        expand("out/kraken2/{sample}/{sample}.tsv",sample=config["samples"]),
        expand("out/kraken2/{sample}/{sample}_bracken_species.txt",sample=config["samples"]),

        #SPADES
        expand("out/scaffolds/{sample}_scaffolds.fasta",sample=config["samples"]),

        #QUAST
        expand("out/quast/{sample}/report.html",sample=config["samples"]),

        #MLST
        expand("out/mlst/{sample}_mlst.tsv",sample=config["samples"]),

        # KLEBORATE
        expand("out/kleborate/{sample}",sample=config["samples"]),

        # COVERAGE
        expand("out/coverage/{sample}_sequencing_coverage.txt",sample=config["samples"]),

        #BAKTA
        expand("out/bakta/{sample}/{sample}.faa", sample=config["samples"]),
        expand("out/bakta/{sample}/{sample}.fna",sample=config["samples"]),
        expand("out/bakta/{sample}/{sample}.gbff",sample=config["samples"]),
        expand("out/bakta/{sample}/{sample}.gff3",sample=config["samples"]),

        #NCBI AMRFINDER PLUS
        expand("out/amrfinder/{sample}.tsv", sample=config["samples"]),

        # PLASMIDFINDER
        expand("out/plasmids/{sample}/results_tab.tsv",sample=config["samples"]),


rule plasmid_finder:
    input:
        fasta="out/prokka/{sample}/{sample}.fna",
    output:
        out_dir=directory("out/plasmids/{sample}"),
        plasmids="out/plasmids/{sample}/results_tab.tsv",
    shell:
        """
        mkdir -p out/plasmids/{wildcards.sample}

        #plasmidfinder.py -i {input.fasta} -o {output.out_dir} -x -p ~/plasmidfinder/plasmidfinder_db

        # I added plasmidfinder to ~/.profile
            python3 ~/plasmidfinder/plasmidfinder.py -i {input.fasta} -o {output.out_dir} -x -p ~/plasmidfinder/plasmidfinder_db
        """

rule amrfinder:
    input:
        fna="out/bakta/{sample}/{sample}.fna",
        faa="out/bakta/{sample}/{sample}.faa",
        gff="out/bakta/{sample}/{sample}.gff3",
        mlst="out/mlst/{sample}_mlst.tsv",
    output:
        amr="out/amrfinder/{sample}.tsv",
    threads: config["threads"]
    shell:
        """
        SAMPLE_STRAIN=$(gawk '{{print $2}}' {input.mlst})
        
        # Set organism-specific flag if available
        ORGANISM_FLAG=""

        case $SAMPLE_STRAIN in
            ecoli*              )   ORGANISM_FLAG="-O Escherichia";;
            klebsiella          )   ORGANISM_FLAG="-O Klebsiella_pneumoniae";;
            paeruginosa         )   ORGANISM_FLAG="-O Pseudomonas_aeruginosa";;
        neisseria           )   ORGANISM_FLAG="-O Neisseria_gonorrhoeae";;
        kaerogenes          )   ORGANISM_FLAG="-O Klebsiella_pneumoniae";;
        efaecalis           )   ORGANISM_FLAG="-O Enterococcus_faecalis";;
        efaecium            )   ORGANISM_FLAG="-O Enterococcus_faecium";;
        spneumoniae         )   ORGANISM_FLAG="-O Streptococcus_pneumoniae";;
        hinfluenzae         )   ORGANISM_FLAG="-O Haemophilus_influenzae";;
        *                   )   ORGANISM_FLAG="";;
        esac

        mkdir -p out/amrfinder

        amrfinder -a bakta --gff {input.gff} -n {input.fna} -p {input.faa} $ORGANISM_FLAG --plus --threads 1 -d /home/raulc/Desktop/Bakta_db/db/amrfinderplus-db/2026-01-21.1/ > {output.amr}
        """

rule bakta:
    input:
        scaffold = "out/scaffolds/{sample}_scaffolds.fasta",
        mlst = "out/mlst/{sample}_mlst.tsv",
        db = "/home/raulc/Desktop/Bakta_db/db" # HARDCODED FOR NOW
    output:
        faa = "out/bakta/{sample}/{sample}.faa",
        fna = "out/bakta/{sample}/{sample}.fna",
        gbff = "out/bakta/{sample}/{sample}.gbff",
        gff3 = "out/bakta/{sample}/{sample}.gff3",
        out_dir=directory("out/bakta/{sample}"),
    conda:"envs/bakta.yaml"
    threads: config["threads"]
    shell:
        """
        #SAMPLE_STRAIN=$(echo {input.mlst} | gawk '{{print $2}}')
        SAMPLE_STRAIN=$(gawk '{{print $2}}' {input.mlst})


        case $SAMPLE_STRAIN in
        paeruginosa         )   GENUS="Pseudomonas"; SPECIES="aeruginosa";;
        spneumoniae         )   GENUS="Streptococcus"; SPECIES="pneumoniae";; 
            ecoli*              )   GENUS="Escherichia"; SPECIES="coli";;
            klebsiella          )   GENUS="Klebsiella"; SPECIES="pneumoniae";;
            kaerogenes          )   GENUS="Klebsiella"; SPECIES="aerogenes";;
        neisseria           )   GENUS="Neisseria"; SPECIES="gonorrhoeae";;
        efaecalis           )   GENUS="Enterococcus"; SPECIES="faecalis";;
        efaecium            )   GENUS="Enterococcus"; SPECIES="faecium";;
        -                   )   GENUS="Legionella"; SPECIES="pneumophila";;
        hinfluenzae         )   GENUS="Haemophilus"; SPECIES="influenzae";;
        esac
        
        mkdir -p out/bakta/{wildcards.sample}
        bakta --db {input.db} --genus $GENUS --species $SPECIES --compliant --prefix {wildcards.sample} --output {output.out_dir} -t 8 -f {input.scaffold}
        """


rule calculate_sequencing_coverage:
    input:
        bam="out/mapping/{sample}.sorted.bam",
        bai="out/mapping/{sample}.sorted.bam.bai",
    output:
        coverage="out/coverage/{sample}_sequencing_coverage.txt",
    shell:
        """
        mkdir -p out/coverage

        AVG_DEPTH=$(samtools depth -a {input.bam} | \
            awk '{{sum+=$3}} END {{print sum/NR}}')

        echo -e "Average sequencing coverage:\t${{AVG_DEPTH}}x" > {output.coverage}
        samtools coverage {input.bam} >> {output.coverage}
        """

rule map_reads_to_assembly:
    input:
        assembly="out/scaffolds/{sample}_scaffolds.fasta",
        r1="all_data/{sample}_R1_001.fastq.gz",
        r2="all_data/{sample}_R2_001.fastq.gz",
    output:
        bam="out/mapping/{sample}.sorted.bam",
        bai="out/mapping/{sample}.sorted.bam.bai",
    threads: 4
    shell:
        """
        mkdir -p out/mapping

        bwa index {input.assembly}
        bwa mem -t {threads} {input.assembly} {input.r1} {input.r2} | \
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule kleborate:
    input:
        scaffold="scaffolds/{sample}_scaffolds.fasta",
    output:
        out_dir=directory("out/kleborate/{sample}"),
    params:
        amr_p="klebsiella_pneumo_complex__amr",
        mlst_p="klebsiella_pneumo_complex__mlst",
        kaptive_p="klebsiella_pneumo_complex__kaptive",
        ybst="klebsiella__ybst",
        cbst="klebsiella__cbst",
        abst="klebsiella__abst",
        smst="klebsiella__smst",
        rmst="klebsiella__rmst",
        rmpa2="klebsiella__rmpa2",
        virulence="klebsiella_pneumo_complex__virulence_score",
    shell:
        """
        mkdir -p out/kleborate/{wildcards.sample}
        kleborate -a {input.scaffold} -m {params.mlst_p} -o {output.out_dir}
        kleborate -a {input.scaffold} -m {params.amr_p} -o {output.out_dir}
        kleborate -a {input.scaffold} -m {params.kaptive_p} -o {output.out_dir}
        kleborate -a {input.scaffold} -m {params.ybst} -o {output.out_dir}
        kleborate -a {input.scaffold} -m {params.cbst} -o {output.out_dir}
        kleborate -a {input.scaffold} -m {params.abst} -o {output.out_dir}
        kleborate -a {input.scaffold} -m {params.smst} -o {output.out_dir}
        kleborate -a {input.scaffold} -m {params.rmst} -o {output.out_dir}
        kleborate -a {input.scaffold} -m {params.rmpa2} -o {output.out_dir}
        kleborate -a {input.scaffold} -m {params.virulence} -o {output.out_dir}
        """

rule mlst:
    input:
        scaffold="out/scaffolds/{sample}_scaffolds.fasta",
    output:
        file="out/mlst/{sample}_mlst.tsv",
    threads: config["threads"]
    shell:
        "mlst --threads {threads} --quiet --label {wildcards.sample}_mlst {input.scaffold} > {output.file}"

rule quast:
    input:
        scaffold="out/scaffolds/{sample}_scaffolds.fasta",
    output:
        html="out/quast/{sample}/report.html",
    shell:
        "mkdir -p out/quast/{wildcards.sample} && "
        "quast.py -o out/quast/{wildcards.sample} {input.scaffold}"

rule copy_scaffolds:
    input:
        scaffolds="out/spades/{sample}/scaffolds.fasta",
    output:
        copied_scaffolds="out/scaffolds/{sample}_scaffolds.fasta",
    shell:
        "cp {input.scaffolds} {output.copied_scaffolds}"


rule spades:
    input:
        read1="out/trim/{sample}_val_1.fq.gz",
        read2="out/trim/{sample}_val_2.fq.gz",
    output:
        out_dir=directory("out/spades/{sample}"),
        scaffolds="out/spades/{sample}/scaffolds.fasta",
    threads: config["threads"]
    shell:
        "spades.py --isolate -1 {input.read1} -2 {input.read2} -o {output.out_dir} -t {threads}"


rule kraken:
    input:
        read1="data_unzipped/{sample}_val_1.fq",
        read2="data_unzipped/{sample}_val_2.fq",
        r1="data_unzipped/{sample}_L001_R1_001.fastq",
        r2="data_unzipped/{sample}_L001_R2_001.fastq",
        kraken_db="/home/raulc/Desktop/work/Sequenziamenti/20251015_corsa/kraken_db",
    output:
        kreport="out/kraken2/{sample}/{sample}.kreport",
        out="out/kraken2/{sample}/{sample}.tsv",
        bracken_out="out/kraken2/{sample}/{sample}_bracken_species.txt",
    threads: config["threads"]
    shell:
        """
        mkdir -p out/kraken2/{wildcards.sample}
        kraken2 --db {input.kraken_db} --paired --use-names {input.read1} {input.read2} --output {output.out} --report {output.kreport} --threads {threads}
        bracken -d {input.kraken_db} -i {output.kreport} -o {output.bracken_out}
        """


rule gunzip_fastq:
    input:
        r1="all_data/{sample}_R1_001.fastq.gz",
        r2="all_data/{sample}_R2_001.fastq.gz",
        r3="out/trim/{sample}_val_1.fq.gz",
        r4="out/trim/{sample}_val_2.fq.gz",
    output:
        r1="data_unzipped/{sample}_L001_R1_001.fastq",
        r2="data_unzipped/{sample}_L001_R2_001.fastq",
        r3="data_unzipped/{sample}_val_1.fq",
        r4="data_unzipped/{sample}_val_2.fq",
    shell:
        """
        mkdir -p data_unzipped
        gunzip -c {input.r1} > {output.r1}
        gunzip -c {input.r2} > {output.r2}
        gunzip -c {input.r3} > {output.r3}
        gunzip -c {input.r4} > {output.r4}
        """


rule fastqc_trimmed:
    input:
        read1="out/trim/{sample}_val_1.fq.gz",
        read2="out/trim/{sample}_val_2.fq.gz",
    output:
        html1="out/fastqc/trimmed/{sample}/{sample}_val_1_fastqc.html",
        html2="out/fastqc/trimmed/{sample}/{sample}_val_2_fastqc.html",
        zip1="out/fastqc/trimmed/{sample}/{sample}_val_1_fastqc.zip",
        zip2="out/fastqc/trimmed/{sample}/{sample}_val_2_fastqc.zip",
    threads: config["threads"]
    shell:
        "mkdir -p out/fastqc/trimmed/{wildcards.sample} && "
        "fastqc -t {threads} -o out/fastqc/trimmed/{wildcards.sample} {input.read1} {input.read2}"


rule trim_galore:
    input:
        read1="all_data/{sample}_R1_001.fastq.gz",
        read2="all_data/{sample}_R2_001.fastq.gz",
    output:
        trimmed_r1="out/trim/{sample}_val_1.fq.gz",
        trimmed_r2="out/trim/{sample}_val_2.fq.gz",
    params:
        quality=config["trimming"]["quality"],
        length=config["trimming"]["length"],
    threads: config["threads"]
    shell:
        "mkdir -p out/trim && "
        "trim_galore -q {params.quality} --length {params.length} --trim-n "
        "--basename {wildcards.sample} -o out/trim --paired {input.read1} {input.read2} --cores {threads}"


rule fastqc_reads:
    input:
        read1="all_data/{sample}_R1_001.fastq.gz",
        read2="all_data/{sample}_R2_001.fastq.gz",
    output:
        html1="out/fastqc/raw/{sample}/{sample}_R1_001_fastqc.html",
        html2="out/fastqc/raw/{sample}/{sample}_R2_001_fastqc.html",
        zip1="out/fastqc/raw/{sample}/{sample}_R1_001_fastqc.zip",
        zip2="out/fastqc/raw/{sample}/{sample}_R2_001_fastqc.zip"
    threads: config["threads"]
    shell:
        "mkdir -p out/fastqc/raw/{wildcards.sample} && "
        "fastqc -t {threads} -o out/fastqc/raw/{wildcards.sample} {input.read1} {input.read2}"

