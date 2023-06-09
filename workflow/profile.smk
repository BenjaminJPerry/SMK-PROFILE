# 2023 Benjamin J Perry
# MIT License
# Copyright (c) 2023 Benjamin J Perry
# Version: 1.0
# Maintainer: Benjamin J Perry
# Email: ben.perry@agresearch.co.nz

configfile: "config/config.yaml"


import os
import pandas as pd


def get_passing_FIDs(seqkitOut):
    import pandas as pd
    qc_stats = pd.read_csv(seqkitOut, delimiter = "\s+")
    qc_stats["num_seqs"] = qc_stats["num_seqs"].str.replace(",", "").astype(int)
    qc_passed = qc_stats.loc[qc_stats["num_seqs"].astype(int) > 50000]
    return qc_passed['file'].str.split("/").str[-1].str.split(".").str[0].tolist()


FIDs = FIDs = get_passing_FIDs("results/00_QC/seqkit.report.KDR.txt")


onstart:
    print(f"Working directory: {os.getcwd()}")
    print("TOOLS: ")
    os.system('echo "  bash: $(which bash)"')
    os.system('echo "  PYTHON: $(which python)"')
    os.system('echo "  CONDA: $(which conda)"')
    os.system('echo "  SNAKEMAKE: $(which snakemake)"')
    print(f"Env TMPDIR = {os.environ.get('TMPDIR', '<n/a>')}")
    os.system('echo "  PYTHON VERSION: $(python --version)"')
    os.system('echo "  CONDA VERSION: $(conda --version)"')


rule all:
    input:
        "results/kraken2.counts.tsv",
        "results/kraken2.counts.biom",
        "results/bracken.k2.counts.tsv",
#       "results/bracken.k2.counts.biom",
        #"results/centrifuge.counts.tsv",
        #"results/centrifuge.counts.biom",
        #expand("results/03_humann3Uniref50EC/{sample}_pathcoverage.tsv", sample=FIDs),


localrules:
    generateCentrifugeSampleSheet,


rule generateCentrifugeSampleSheet:
    output:
        sampleSheet = "resources/centrifugeSampleSheet.tsv",
    threads: 2
    shell:
        "./workflow/scripts/generate_centrifuge_sample_sheet.sh -d results/02_kneaddata -p fastq.gz -o {output.sampleSheet} "


wildcard_constraints:
    sample="\w+"


rule centrifugeGTDB:
    input:
        sampleSheet = "resources/centrifugeSampleSheet.tsv",
#        KDRs = "results/02_kneaddata/{sample}.fastq",
    output:
#        out = "results/03_centrifuge/{sample}.GTDB.centrifuge",
#        report = "results/03_centrifuge/{sample}.GTDB.centrifuge.report",
        out = expand("results/03_centrifuge/{sample}.GTDB.centrifuge", sample = FIDs),
        report = expand("results/03_centrifuge/{sample}.GTDB.centrifuge.report", sample = FIDs),
    log:
        "logs/centrifugeGTDB.log",
#        "logs/centrifugeGTDB.{sample}.log",
    benchmark:
        "benchmarks/centrifugeGTDB.txt"
#        "benchmarks/centrifugeGTDB.{sample}.txt"
    conda:
        "centrifuge"
    threads: 64
    resources:
        mem_gb = lambda wildcards, attempt: 160 + ((attempt - 1) * 20),
        time = "05-00:00:00",
        partition = "milan"
    shell:
        "centrifuge "
        "-x /nesi/nobackup/agresearch03843/centrifuge/centrifuge/GTDB "
        "--sample-sheet {input.sampleSheet} "
#        "-U {input.KDRs} "
#        "--report {output.report} "
#        "-S {output.out} "
        "-t "
        "--threads {threads} "
        "2>&1 | tee {log}"


rule centrifugeKrakenReport:
    input:
        centrifuge = "results/03_centrifuge/{sample}.GTDB.centrifuge",
    output:
        centrifugeKraken2 = "results/03_centrifuge/{sample}.centrifuge",
    log:
        "logs/centrifugeKrakenReport/{sample}.centrifuge.to.kraken2.log",
    benchmark:
        "benchmarks/centrifugeKrakenReport.{sample}.txt"
    conda:
        "centrifuge"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 2 + ((attempt - 1) + 2),
        time = lambda wildcards, attempt: 2 + ((attempt - 1) + 2),
        partition = "large,milan"
    shell:
        "centrifuge-kreport "
        "-x /nesi/nobackup/agresearch03843/centrifuge/centrifuge/GTDB "
        "{input.centrifuge} > "
        "{output.centrifugeKraken2}"


rule taxpastaCentrifugeTable:
    input:
        expand("results/03_centrifuge/{sample}.GTDB.centrifuge", sample = FIDs),
    output:
        "results/centrifuge.counts.tsv",
    benchmark:
        "benchmarks/taxpastaCentrifugeTable.txt"
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "large,milan"
    shell:
        "taxpasta merge "
        "-p centrifuge "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /nesi/nobackup/agresearch03843/centrifuge/centrifuge "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at genus "
        "{input} "


rule taxpastaCentrifugeBiom:
    input:
        expand("results/03_centrifuge/{sample}.GTDB.centrifuge", sample = FIDs),
    output:
        "results/centrifuge.counts.biom",
    benchmark:
        "benchmarks/taxpastaCentrifugeBiom.txt"
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "large,milan"
    shell:
        "taxpasta merge "
        "-p centrifuge "
        "-o {output} "
        "--output-format BIOM "
        "--taxonomy /nesi/nobackup/agresearch03843/centrifuge/centrifuge "
        "--add-name "
        "--summarise-at genus "
        "{input} "


#KRAKEN2 RULES
rule kraken2GTDB:
    input:
        KDRs = "results/02_kneaddata/{sample}.fastq.gz",
    output:
        k2OutputGTDB = "results/03_kraken2GTDB/{sample}.k2",
        k2ReportGTDB = "results/03_kraken2GTDB/{sample}.kraken2",
    log:
        "logs/kraken2GTDB.{sample}.GTDB.log",
    benchmark:
        "benchmarks/kraken2GTDB.{sample}.txt"
    conda:
        "kraken2"
    threads: 32
    resources:
        mem_gb = lambda wildcards, attempt: 340 + ((attempt - 1) * 20),
        time = lambda wildcards, attempt: 15 + ((attempt - 1) * 60),
        partition = "milan"
    shell:
        "kraken2 "
        "--use-names "
        # "--quick "
        "--db /nesi/nobackup/agresearch03843/kraken2-GTDB/GTDB "
        "-t {threads} "
        "--report {output.k2ReportGTDB} "
        "--report-minimizer-data "
        "{input.KDRs} "
        "--output {output.k2OutputGTDB} "
        "2>&1 | tee {log} "


rule taxpastaKraken2:
    input:
        expand("results/03_kraken2GTDB/{sample}.kraken2", sample = FIDs),
    output:
        "results/metagenomes.kraken2.genus.GTDB207.tsv",
    benchmark:
        "benchmarks/taxpastaKraken2.txt"
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "large,milan"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format TSV "
        "--taxonomy /nesi/nobackup/agresearch03843/kraken2-GTDB/GTDB/taxonomy "
        "--add-name "
        "--add-rank "
        "--add-lineage "
        "--summarise-at genus "
        "{input} "


rule taxpastaKraken2Biom:
    input:
        expand("results/03_kraken2GTDB/{sample}.kraken2", sample = FIDs),
    output:
        "results/kraken2.counts.biom",
    benchmark:
        "benchmarks/taxpastaKraken2Biom.txt"
    conda:
        "taxpasta"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 1440 + ((attempt - 1) * 1440),
        partition = "large,milan"
    shell:
        "taxpasta merge "
        "-p kraken2 "
        "-o {output} "
        "--output-format BIOM "
        "--taxonomy /nesi/nobackup/agresearch03843/kraken2-GTDB/GTDB/taxonomy  " 
        "--add-name "
        "--summarise-at genus "
        "{input} "


# BRACKEN RULES
rule brackenSpecies:
    input:
        k2ReportGTDB = "results/03_kraken2GTDB/{sample}.kraken2",
    output:
        bOutput = "results/03_brackenSpecies/{sample}.bracken",
        bReport = "results/03_brackenSpecies/{sample}.br",
    log:
        "logs/brackenSpecies.{sample}.log",
    benchmark:
        "benchmarks/brackenSpecies.{sample}.txt"
    conda:
        "kraken2"
    threads: 2
    resources:
        mem_gb = 1,
        time = 2,
        partition = "large,milan"
    shell:
        "bracken "
        "-d /nesi/nobackup/agresearch03843/kraken2-GTDB/GTDB "
        "-i {input.k2ReportGTDB} "
        "-o {output.bOutput} "
        "-w {output.bReport} "
        "-r 80 "
        "-l S "
        "-t 1 "
        "&> {log} "


rule mergeKraken2Bracken:
    input:
        expand("results/03_brackenSpecies/{sample}.bracken", sample = FIDs),
    output:
        "results/metagenomes.bracken.k2.species.GTDB207.tsv",
    benchmark:
        "benchmarks/taxpastaKraken2Bracken.txt"
    conda:
        "kraken2"
    threads: 2
    resources:
        mem_gb = lambda wildcards, attempt: 48 + ((attempt - 1) * 8),
        time = lambda wildcards, attempt: 2880 + ((attempt - 1) * 1440),
        partition = "large,milan"
    shell:
        "combine_bracken_outputs.py "
        "--files {input} "
        "-o {output}"


# HUMANN# RULES
rule humann3Uniref50EC:
    input:
        kneaddataReads = "results/02_kneaddata/{sample}.fastq.gz",
    output:
        genes = "results/03_humann3Uniref50EC/{sample}_genefamilies.tsv",
        pathways = "results/03_humann3Uniref50EC/{sample}_pathabundance.tsv",
        pathwaysCoverage = "results/03_humann3Uniref50EC/{sample}_pathcoverage.tsv",
    log:
        "logs/humann3.{sample}.uniref50EC.log",
    benchmark:
        "benchmarks/humann3Uniref50EC.{sample}.txt"
    conda:
        "biobakery"
    threads: 24
    resources:
        mem_gb = lambda wildcards, attempt: 32 + ((attempt - 1) + 8),
        time = lambda wildcards, attempt: 30 + ((attempt - 1) + 30),
        partition = "large,milan"
    message:
        "humann3 profiling with uniref50EC: {wildcards.sample}\n"
    shell:
        "humann3 "
        "--memory-use maximum "
        "--threads {threads} "
        "--bypass-nucleotide-search "
        "--search-mode uniref50 "
        "--protein-database /nesi/nobackup/agresearch03843/biobakery/biobakery/humann3/unirefECFilt "
        "--input-format fastq.gz "
        "--output results/03_humann3Uniref50EC "
        "--input {input.kneaddataReads} "
        "--output-basename {wildcards.sample} "
        "--o-log {log} "
        "--remove-temp-output "
