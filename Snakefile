import os
from glob import glob

configfile: "config/config.yaml"

# dataset is passed via --config dataset=dataset1
DATASET = config["dataset"]

DATASET_DIR = f"data/{DATASET}"
RESULTS_DIR = f"results/{DATASET}"

short = [l.strip() for l in open(DATASET_DIR+"/short_reads.txt") if l.strip()]
long  = [l.strip() for l in open(DATASET_DIR+"/long_reads.txt") if l.strip()]

assert len(short) == len(long), "Short/long accession count mismatch"

SAMPLES = [f"S{i:02d}" for i in range(len(short))]

SHORT_ACC = dict(zip(SAMPLES, short))
LONG_ACC  = dict(zip(SAMPLES, long))

############################################
# Rule order / final targets
############################################

rule all:
    input:
        # short-read-only assemblies
        expand(
            f"{RESULTS_DIR}/assemblies/short/{{sample}}/contigs.fasta",
            sample=SAMPLES
        ),
        # long-read-only assemblies
        expand(
            f"{RESULTS_DIR}/assemblies/long/{{sample}}/contigs.fasta",
            sample=SAMPLES
        ),
        # hybrid assemblies
        expand(
            f"{RESULTS_DIR}/assemblies/hybrid/{{sample}}/contigs.fasta",
            sample=SAMPLES
        ),
        # multi assembly
        f"{RESULTS_DIR}/assemblies/multi/contigs.fasta",

        # binning
        expand(f"{RESULTS_DIR}/bins/coassembly/{{tool}}/bins.done",tool=config["binning_tools"]),
        expand(f"{RESULTS_DIR}/bins/single/{{tool}}/{{assembly_type}}/{{sample}}/bins.done", tool=config["binning_tools"], assembly_type=["short", "long", "hybrid"], sample=SAMPLES),
        expand(f"{RESULTS_DIR}/bins/multi/{{tool}}/{{assembly_type}}/bins.done", tool=config["binning_tools"], assembly_type=["short", "long", "hybrid"]),
        # QC
        # Single-sample QC
        expand(f"{RESULTS_DIR}/qc/single/{{assembly_type}}/{{sample}}/qc.done",
               assembly_type=["short","long","hybrid"], sample=SAMPLES),
        # Multi-sample QC
        f"{RESULTS_DIR}/qc/multi/short/qc.done",

############################################
# 1. Download SRA / dump files
############################################

rule download_short_reads:
    output:
        sra=f"{RESULTS_DIR}/raw/short/{{sample}}.sra"
    params:
        acc=lambda wc: SHORT_ACC[wc.sample]
    log:
        f"logs/{DATASET}/download/short/{{sample}}.log"
    shell:
        """
        prefetch {params.acc} --output-file {output.sra}
        """

rule download_long_reads:
    output:
        sra=f"{RESULTS_DIR}/raw/long/{{sample}}.sra"
    params:
        acc=lambda wc: LONG_ACC[wc.sample]
    log:
        f"logs/{DATASET}/download/long/{{sample}}.log"
    shell:
        """
        prefetch {params.acc} --output-file {output.sra}
        """

############################################
# 2. Extract FASTQ
############################################

rule sra_to_fastq_short:
    input:
        sra=f"{RESULTS_DIR}/raw/short/{{sample}}.sra"
    output:
        r1=f"{RESULTS_DIR}/fastq/short/{{sample}}_R1.fq.gz",
        r2=f"{RESULTS_DIR}/fastq/short/{{sample}}_R2.fq.gz"
    params:
        acc=lambda wc: SHORT_ACC[wc.sample]
    log:
        f"logs/{DATASET}/fastq/short/{{sample}}.log"
    shell:
        """
        fasterq-dump {input.sra} --split-files -O {RESULTS_DIR}/fastq/short \
        2> {log}
        gzip -f {RESULTS_DIR}/fastq/short/{wildcards.sample}_1.fastq
        gzip -f {RESULTS_DIR}/fastq/short/{wildcards.sample}_2.fastq
        mv {RESULTS_DIR}/fastq/short/{wildcards.sample}_1.fastq.gz {output.r1}
        mv {RESULTS_DIR}/fastq/short/{wildcards.sample}_2.fastq.gz {output.r2}
        """

rule sra_to_fastq_long:
    input:
        sra=f"{RESULTS_DIR}/raw/long/{{sample}}.sra"
    output:
        fq=f"{RESULTS_DIR}/fastq/long/{{sample}}.fq.gz"
    params:
        acc=lambda wc: LONG_ACC[wc.sample]
    log:
        f"logs/{DATASET}/fastq/long/{{sample}}.log"
    shell:
        """
        fasterq-dump {input.sra} -O {RESULTS_DIR}/fastq/long 2> {log}
        gzip -f {RESULTS_DIR}/fastq/long/{wildcards.sample}.fastq
        mv {RESULTS_DIR}/fastq/long/{wildcards.sample}.fastq.gz {output.fq}
        """

############################################
# 3. Preprocessing
############################################

rule preprocess_short:
    input:
        r1=f"{RESULTS_DIR}/fastq/short/{{sample}}_R1.fq.gz",
        r2=f"{RESULTS_DIR}/fastq/short/{{sample}}_R2.fq.gz"
    output:
        r1=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.fq.gz",
        r2=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.fq.gz"
    log:
        f"logs/{DATASET}/preprocess/short/{{sample}}.log"
    shell:
        """
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          > {log} 2>&1
        """

rule preprocess_long:
    input:
        fq=f"{RESULTS_DIR}/fastq/long/{{sample}}.fq.gz"
    output:
        fq=f"{RESULTS_DIR}/preprocess/long/{{sample}}.fq.gz"
    log:
        f"logs/{DATASET}/preprocess/long/{{sample}}.log"
    shell:
        """
        # For long reads, fastp can still be used, but only single input
        fastp \
          -i {input.fq} -o {output.fq} \
          > {log} 2>&1
        """

############################################
# 4. Assemblies
############################################

# 4.1 Single-sample assemblies

rule assemble_single_short:
    input:
        r1=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.fq.gz",
        r2=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.fq.gz"
    output:
        f"{RESULTS_DIR}/assemblies/single/short/{{sample}}/contigs.fasta"
    threads: config["threads"]
    log:
        f"logs/{DATASET}/assembly/single/short/{{sample}}.log"
    shell:
        """
        spades.py \
          -1 {input.r1} -2 {input.r2} \
          -o {RESULTS_DIR}/assemblies/single/short/{wildcards.sample} \
          > {log} 2>&1
        """

rule assemble_single_long:
    input:
        fq=f"{RESULTS_DIR}/preprocess/long/{{sample}}.fq.gz"
    output:
        f"{RESULTS_DIR}/assemblies/single/long/{{sample}}/contigs.fasta"
    threads: config["threads"]
    log:
        f"logs/{DATASET}/assembly/single/long/{{sample}}.log"
    shell:
        """
        spades.py \
          --pacbio {input.fq} \
          -o {RESULTS_DIR}/assemblies/single/long/{wildcards.sample} \
          > {log} 2>&1
        """

rule assemble_hybrid:
    input:
        r1=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.fq.gz",
        r2=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.fq.gz",
        long=f"{RESULTS_DIR}/preprocess/long/{{sample}}.fq.gz"
    output:
        f"{RESULTS_DIR}/assemblies/hybrid/{{sample}}/contigs.fasta"
    threads: config["threads"]
    log:
        f"logs/{DATASET}/assembly/hybrid/{{sample}}.log"
    shell:
        """
        spades.py \
          -1 {input.r1} -2 {input.r2} \
          --pacbio {input.long} \
          -o {RESULTS_DIR}/assemblies/hybrid/{wildcards.sample} \
          > {log} 2>&1
        """

# 4.2 Multi-sample assemblies (dataset-level)

rule assemble_multi_short:
    input:
        r1=expand(f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.fq.gz", sample=SAMPLES),
        r2=expand(f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.fq.gz", sample=SAMPLES)
    output:
        f"{RESULTS_DIR}/assemblies/multi/short/contigs.fasta"
    threads: config["threads"]
    log:
        f"logs/{DATASET}/assembly/multi_short.log"
    shell:
        """
        megahit \
          -1 {','.join(input.r1)} \
          -2 {','.join(input.r2)} \
          -o {RESULTS_DIR}/assemblies/multi/short \
          > {log} 2>&1
        """

############################################
# 5. Binning
############################################

rule bin_coassembly:
    input:
        assembly=f"{RESULTS_DIR}/assemblies/multi/short/contigs.fasta"
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/{{tool}}/bins.done")
    params:
        tool=lambda wc: wc.tool
    threads: config["threads"]
    log:
        f"logs/{DATASET}/binning/coassembly/{{tool}}.log"
    shell:
        """
        # example placeholder, replace with actual binning command
        binning_tool --input {input.assembly} --output {wildcards.tool}_bins \
        > {log} 2>&1
        touch {output}
        """

rule bin_single_sample:
    input:
        assembly=lambda wc: f"{RESULTS_DIR}/assemblies/{wc.assembly_type}/{wc.sample}/contigs.fasta"
    output:
        touch(f"{RESULTS_DIR}/bins/{{tool}}/{{assembly_type}}/{{sample}}/bins.done")
    params:
        tool=lambda wc: wc.tool
    threads: config["threads"]
    log:
        f"logs/{DATASET}/binning/single/{{tool}}/{{assembly_type}}/{{sample}}.log"
    shell:
        """
        binning_tool --input {input} --output {wildcards.tool}_{wildcards.sample}_{wildcards.assembly_type}_bins \
        > {log} 2>&1
        touch {output}
        """

rule bin_multi_sample:
    input:
        assemblies=lambda wc: expand(
            f"{RESULTS_DIR}/assemblies/single/{wc.assembly_type}/{{sample}}/contigs.fasta",
            sample=SAMPLES
        )
    output:
        touch(f"{RESULTS_DIR}/bins/multi/{{tool}}/{{assembly_type}}/bins.done")
    params:
        tool=lambda wc: wc.tool
    threads: config["threads"]
    log:
        f"logs/{DATASET}/binning/multi/{{tool}}/{{assembly_type}}.log"
    shell:
        """
        binning_tool --input {','.join(input.assemblies)} --output {wildcards.tool}_{wildcards.assembly_type}_multi_bins \
        > {log} 2>&1
        touch {output}
        """

############################################
# 6. Quality control
############################################

# Single-sample QC (short, long, hybrid)
rule qc_single:
    input:
        assembly=lambda wc: f"{RESULTS_DIR}/assemblies/single/{wc.assembly_type}/{wc.sample}/contigs.fasta"
    output:
        touch(f"{RESULTS_DIR}/qc/single/{{assembly_type}}/{{sample}}/qc.done")
    log:
        f"logs/{DATASET}/qc/single/{{assembly_type}}/{{sample}}.log"
    shell:
        """
        quast {input.assembly} -o qc_tmp_{wildcards.sample} > {log} 2>&1
        # you could move or copy results if needed
        touch {output}
        """

# Multi-sample QC (multi_short)
rule qc_multi:
    input:
        assembly=f"{RESULTS_DIR}/assemblies/multi/short/contigs.fasta"
    output:
        touch(f"{RESULTS_DIR}/qc/multi/short/qc.done")
    log:
        f"logs/{DATASET}/qc/multi_short.log"
    shell:
        """
        quast {input.assembly} -o qc_tmp_multi_short > {log} 2>&1
        touch {output}
        """

