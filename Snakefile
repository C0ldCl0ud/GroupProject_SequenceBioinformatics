wildcard_constraints:
    sample = "S\\d+"

import os
from glob import glob

configfile: "config/config.yaml"

# dataset is passed via --config dataset=dataset1
DATASET = config["dataset"]
DATASET_DIR = f"data/{DATASET}"
RESULTS_DIR = f"results/{DATASET}"

short = [l.strip() for l in open(DATASET_DIR+"/short_reads.txt") if l.strip()]
long  = [l.strip() for l in open(DATASET_DIR+"/long_reads.txt") if l.strip()]

SAMPLES = [f"S{i:02d}" for i in range(len(short))]

#check long read technology
LONG_TECH = config["long_read_technology"].lower()
if LONG_TECH == "nanopore":
    LONG_PREPROCESSED = f"{RESULTS_DIR}/preprocess/long/{{sample}}.final.fq.gz"
elif LONG_TECH == "pacbio":
    LONG_PREPROCESSED = f"{RESULTS_DIR}/fastq/long/{{sample}}.fq.gz"
else:
    raise ValueError("Unsupported long_read_technology")

#check human gut removal
if bool(config.get("remove_human_reads", False)):
    SHORT_FINAL_R1 = f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.nohuman.fq.gz"
    SHORT_FINAL_R2 = f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.nohuman.fq.gz"
    LONG_FINAL = f"{RESULTS_DIR}/preprocess/long/{{sample}}.nohuman.fq.gz"
else:
    SHORT_FINAL_R1 = f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.fastp.fq.gz"
    SHORT_FINAL_R2 = f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.fastp.fq.gz"
    LONG_FINAL = LONG_PREPROCESSED

#assembly modes
if LONG_TECH == "pacbio":
    FLYE_MODE = "--pacbio-hifi"
elif LONG_TECH == "nanopore":
    FLYE_MODE = "--nano-raw"

#assemply post processing with pilon
if LONG_TECH == "nanopore":
    HYBRID_FINAL = f"{RESULTS_DIR}/assemblies/single/hybrid/{{sample}}/contigs.pilon.fasta"
else:
    HYBRID_FINAL = f"{RESULTS_DIR}/assemblies/single/hybrid/{{sample}}/contigs.fasta"



assert len(short) == len(long), "Short/long accession count mismatch"

SHORT_ACC = dict(zip(SAMPLES, short))
LONG_ACC  = dict(zip(SAMPLES, long))

############################################
# Rule order / final targets
############################################

rule all:
    input:
        [
        #expand(f"{RESULTS_DIR}/raw/short/{{sample}}", sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/raw/long/{{sample}}", sample=SAMPLES),

        #expand(f"{RESULTS_DIR}/fastq/short/{{sample}}_R1.fq.gz", sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/fastq/short/{{sample}}_R2.fq.gz", sample=SAMPLES),
            #expand(f"{RESULTS_DIR}/fastq/long/{{sample}}.fq.gz", sample=SAMPLES),

            #expand(f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.fastp.fq.gz", sample=SAMPLES),
            #expand(f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.fastp.fq.gz", sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/preprocess/long/{{sample}}.qcat.fq.gz", sample=SAMPLES),

        #expand(f"{RESULTS_DIR}/preprocess/long/{{sample}}.filtlong1.fq.gz", sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/preprocess/long/{{sample}}.porechop.fq.gz", sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/preprocess/long/{{sample}}.final.fq.gz", sample=SAMPLES),
        
        #expand(f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.nohuman.fq.gz", sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.nohuman.fq.gz", sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/preprocess/long/{{sample}}.nohuman.fq.gz", sample=SAMPLES),

        # single-sample assemblies
        #expand(
        #    f"{RESULTS_DIR}/assemblies/single/{{asm_type}}/{{sample}}/assembly.fasta",
        #    #asm_type=["short","long"],
        #    asm_type=["short","long","hybrid"],
        #    sample=SAMPLES
        #),

        # multi-sample assembly
        #f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta"


        # binning
        expand(f"{RESULTS_DIR}/bins/coassembly/{{tool}}/bins.done",tool=config["binning_tools"]),
        expand(f"{RESULTS_DIR}/bins/single/{{tool}}/{{assembly_type}}/{{sample}}/bins.done", tool=config["binning_tools"], assembly_type=["short", "long", "hybrid"], sample=SAMPLES),
        expand(f"{RESULTS_DIR}/bins/multi/{{tool}}/{{assembly_type}}/{{sample}}/bins.done", tool=config["binning_tools"], sample=SAMPLES, assembly_type=["short", "long", "hybrid"])
        # QC
        # Single-sample QC
        #expand(f"{RESULTS_DIR}/qc/single/{{assembly_type}}/{{sample}}/qc.done",
               #assembly_type=["short","long","hybrid"], sample=SAMPLES),
        # Multi-sample QC
        #f"{RESULTS_DIR}/qc/multi/short/qc.done",
        ]

############################################
# 1. Download SRA / dump files
############################################

rule download_short_reads:
    output:
        sra=directory(f"{RESULTS_DIR}/raw/short/{{sample}}")
    params:
        acc=lambda wc: SHORT_ACC[wc.sample]
    conda:
        "envs/download.yaml"
    log:
        f"logs/{DATASET}/download/short/{{sample}}.log"
    shell:
        """
        prefetch {params.acc} -O {output.sra} > {log} 2>&1
        """

rule download_long_reads:
    output:
        sra_dir = directory(f"{RESULTS_DIR}/raw/long/{{sample}}")
    params:
        acc = lambda wc: LONG_ACC[wc.sample]
    conda:
        "envs/download.yaml"
    log:
        f"logs/{DATASET}/download/long/{{sample}}.log"
    shell:
        """
        prefetch {params.acc} -O {output.sra_dir} > {log} 2>&1
        """

############################################
# 2. Extract FASTQ
############################################
rule sra_to_fastq_short:
    input:
        sra_dir = f"{RESULTS_DIR}/raw/short/{{sample}}"
    output:
        r1 = f"{RESULTS_DIR}/fastq/short/{{sample}}_R1.fq.gz",
        r2 = f"{RESULTS_DIR}/fastq/short/{{sample}}_R2.fq.gz"
    params:
        acc = lambda wc: SHORT_ACC[wc.sample],
        debug = config.get("debug", False),
        n = config.get("debug_reads_short", 200000),
        seed = config.get("debug_seed", 42)
    conda:
        "envs/download.yaml"
    log:
        f"logs/{DATASET}/fastq/short/{{sample}}.log"
    shell:
        """
        fasterq-dump {input.sra_dir}/{params.acc}/{params.acc}.sra \
            --split-files \
            -O {RESULTS_DIR}/fastq/short \
            > {log} 2>&1

        if [ "{params.debug}" = "True" ]; then
            seqtk sample -s{params.seed} \
                {RESULTS_DIR}/fastq/short/{params.acc}_1.fastq {params.n} \
                | gzip -c > {output.r1}

            seqtk sample -s{params.seed} \
                {RESULTS_DIR}/fastq/short/{params.acc}_2.fastq {params.n} \
                | gzip -c > {output.r2}

            rm {RESULTS_DIR}/fastq/short/{params.acc}_1.fastq
            rm {RESULTS_DIR}/fastq/short/{params.acc}_2.fastq
        else
            gzip -f {RESULTS_DIR}/fastq/short/{params.acc}_1.fastq
            gzip -f {RESULTS_DIR}/fastq/short/{params.acc}_2.fastq

            mv {RESULTS_DIR}/fastq/short/{params.acc}_1.fastq.gz {output.r1}
            mv {RESULTS_DIR}/fastq/short/{params.acc}_2.fastq.gz {output.r2}
        fi
        """

rule sra_to_fastq_long:
    input:
        sra_dir = f"{RESULTS_DIR}/raw/long/{{sample}}"
    output:
        fq = f"{RESULTS_DIR}/fastq/long/{{sample}}.fq.gz"
    params:
        acc = lambda wc: LONG_ACC[wc.sample],
        debug = config.get("debug", False),
        n = config.get("debug_reads_long", 5000),
        seed = config.get("debug_seed", 42)
    conda:
        "envs/download.yaml"
    log:
        f"logs/{DATASET}/fastq/long/{{sample}}.log"
    shell:
        """
        fasterq-dump {input.sra_dir}/{params.acc}/{params.acc}.sra \
            -O {RESULTS_DIR}/fastq/long \
            > {log} 2>&1

        if [ "{params.debug}" = "True" ]; then
            seqtk sample -s{params.seed} \
                {RESULTS_DIR}/fastq/long/{params.acc}.fastq {params.n} \
                | gzip -c > {output.fq}

            rm {RESULTS_DIR}/fastq/long/{params.acc}.fastq
        else
            gzip -f {RESULTS_DIR}/fastq/long/{params.acc}.fastq
            mv {RESULTS_DIR}/fastq/long/{params.acc}.fastq.gz {output.fq}
        fi
        """

############################################
# 3. Preprocessing
############################################

rule preprocess_short_fastp:
    input:
        r1=f"{RESULTS_DIR}/fastq/short/{{sample}}_R1.fq.gz",
        r2=f"{RESULTS_DIR}/fastq/short/{{sample}}_R2.fq.gz"
    output:
        r1=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.fastp.fq.gz",
        r2=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.fastp.fq.gz"
    log:
        f"logs/{DATASET}/preprocess/short/{{sample}}.fastp.log"
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          -q 20 \
          --length_required 100 \
          --low_complexity_filter \
          > {log} 2>&1
        """

rule nanopore_qcat:
    input:
        fq=f"{RESULTS_DIR}/fastq/long/{{sample}}.fq.gz"
    output:
        fq=f"{RESULTS_DIR}/preprocess/long/{{sample}}.qcat.fq.gz"
    log:
        f"logs/{DATASET}/preprocess/long/{{sample}}.qcat.log"
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        qcat \
          --trim \
          --detect-middle \
          -i {input.fq} \
          -o {output.fq} \
          > {log} 2>&1
        """

rule nanopore_filtlong_1:
    input:
        fq=f"{RESULTS_DIR}/preprocess/long/{{sample}}.qcat.fq.gz"
    output:
        fq=f"{RESULTS_DIR}/preprocess/long/{{sample}}.filtlong1.fq.gz"
    log:
        f"logs/{DATASET}/preprocess/long/{{sample}}.filtlong1.log"
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        filtlong \
          --min_length 4000 \
          --min_mean_q 80 \
          {input.fq} | gzip > {output.fq} \
          2> {log}
        """

rule nanopore_porechop:
    input:
        fq = f"{RESULTS_DIR}/preprocess/long/{{sample}}.filtlong1.fq.gz"

    output:
        fq=f"{RESULTS_DIR}/preprocess/long/{{sample}}.porechop.fq.gz"
    log:
        f"logs/{DATASET}/preprocess/long/{{sample}}.porechop.log"
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        porechop \
          -i {input.fq} \
          -o {output.fq} \
          --min_split_read_size 4000 \
          > {log} 2>&1
        """

rule nanopore_filtlong_2:
    input:
        fq=f"{RESULTS_DIR}/preprocess/long/{{sample}}.porechop.fq.gz"
    output:
        fq=f"{RESULTS_DIR}/preprocess/long/{{sample}}.final.fq.gz"
    log:
        f"logs/{DATASET}/preprocess/long/{{sample}}.filtlong2.log"
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        filtlong \
          --min_length 4000 \
          --min_mean_q 80 \
          {input.fq} | gzip > {output.fq} \
          2> {log}
        """

rule remove_human_short:
    input:
        r1=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.fastp.fq.gz",
        r2=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.fastp.fq.gz"
    output:
        r1=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.nohuman.fq.gz",
        r2=f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.nohuman.fq.gz"
    threads: config["threads"]
    log:
        f"logs/{DATASET}/preprocess/short/{{sample}}.bowtie2.log"
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        bowtie2 \
          -x {config[hg38_index]} \
          -1 {input.r1} -2 {input.r2} \
          -p {threads} \
          --un-conc-gz {RESULTS_DIR}/preprocess/short/{wildcards.sample}.nohuman.fq.gz \
          > {log} 2>&1

        mv {RESULTS_DIR}/preprocess/short/{wildcards.sample}.nohuman.1.fq.gz {output.r1}
        mv {RESULTS_DIR}/preprocess/short/{wildcards.sample}.nohuman.2.fq.gz {output.r2}
        """

rule remove_human_long:
    input:
        fq=LONG_PREPROCESSED
    output:
        fq=f"{RESULTS_DIR}/preprocess/long/{{sample}}.nohuman.fq.gz"
    threads: config["threads"]
    log:
        f"logs/{DATASET}/preprocess/long/{{sample}}.bowtie2.log"
    conda:
        "envs/preprocessing.yaml"
    shell:
        """
        bowtie2 \
          -x {config[hg38_index]} \
          -U {input.fq} \
          -p {threads} \
          --un-gz {output.fq} \
          > {log} 2>&1
        """

############################################
# 4. Assemblies
############################################

rule assemble_single_short:
    input:
        r1=SHORT_FINAL_R1,
        r2=SHORT_FINAL_R2
    output:
        contigs = f"{RESULTS_DIR}/assemblies/single/short/{{sample}}/assembly.fasta"
    threads: config["threads"]
    log:
        f"logs/{DATASET}/assembly/single/short/{{sample}}.megahit.log"
    conda:
        "envs/assembly.yaml"
    params:
        tmpdir = lambda wc: f"{RESULTS_DIR}/assemblies/single/short/{wc.sample}_tmp"
    shell:
        """
        rm -rf {params.tmpdir}

        megahit \
          -1 {input.r1} \
          -2 {input.r2} \
          -o {params.tmpdir} \
          --min-contig-len 1000 \
          -t {threads} \
          > {log} 2>&1

        mkdir -p $(dirname {output.contigs})
        cp {params.tmpdir}/final.contigs.fa {output.contigs}
        """

rule assemble_single_long:
    input:
        fq = LONG_FINAL
    output:
        contigs = f"{RESULTS_DIR}/assemblies/single/long/{{sample}}/assembly.fasta"
    threads: config["threads"]
    log:
        f"logs/{DATASET}/assembly/single/long/{{sample}}.flye.log"
    conda:
        "envs/assembly.yaml"
    shell:
        """
        flye \
          {FLYE_MODE} {input.fq} \
          --out-dir {RESULTS_DIR}/assemblies/single/long/{wildcards.sample} \
          --threads {threads} \
          --min-overlap 1000 \
          --meta \
          > {log} 2>&1
        """

rule assemble_single_hybrid:
    input:
        r1=SHORT_FINAL_R1,
        r2=SHORT_FINAL_R2,
        contigs=f"{RESULTS_DIR}/assemblies/single/short/{{sample}}/assembly.fasta",
        long=LONG_FINAL
    output:
        contigs=f"{RESULTS_DIR}/assemblies/single/hybrid/{{sample}}/assembly.fasta"
    threads: 16
    log:
        f"logs/{DATASET}/assembly/single/hybrid/{{sample}}.operams.log"
    container:
        #"containers/operams.simg" # remote server
        "containers/operams.sif"   # local machine
    shell:
        """
        set -euo pipefail

        export TMPDIR={RESULTS_DIR}/tmp/{wildcards.sample}
        mkdir -p "$TMPDIR"

        long_unzipped=$(mktemp --suffix=.fastq -p "$TMPDIR")
        gunzip -c {input.long} > "$long_unzipped"

        rm -rf {RESULTS_DIR}/assemblies/single/hybrid/{wildcards.sample}

        perl /operams/OPERA-MS.pl \
          --contig-file {input.contigs} \
          --short-read1 {input.r1} \
          --short-read2 {input.r2} \
          --long-read "$long_unzipped" \
          --no-polishing \
          --out-dir {RESULTS_DIR}/assemblies/single/hybrid/{wildcards.sample} \
          --num-processors {threads} \
          --no-ref-clustering \
          > {log} 2>&1

        cp {RESULTS_DIR}/assemblies/single/hybrid/{wildcards.sample}/contigs.fasta {output.contigs}
        rm -f "$long_unzipped"
        """

rule assemble_coassembly_short:
    input:
        r1=expand(SHORT_FINAL_R1, sample=SAMPLES),
        r2=expand(SHORT_FINAL_R2, sample=SAMPLES)
    output:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta"
    threads: config["threads"]
    log:
        f"logs/{DATASET}/assembly/coassembly/short.megahit.log"
    conda:
        "envs/assembly.yaml"
    params:
        tmpdir = f"{RESULTS_DIR}/assemblies/coassembly/short_tmp",
        r1 = lambda wc, input: ",".join(input.r1),
        r2 = lambda wc, input: ",".join(input.r2)
    shell:
        """
        rm -rf {params.tmpdir}

        megahit \
          -1 {params.r1} \
          -2 {params.r2} \
          -o {params.tmpdir} \
          --min-contig-len 1000 \
          -t {threads} \
          > {log} 2>&1

        mkdir -p $(dirname {output.contigs})
        cp {params.tmpdir}/final.contigs.fa {output.contigs}
        """

############################################
# 5. Binning (mapping + depth + binning)
############################################

############################################
# 5.1 Build indices (ONCE per assembly)
############################################

rule index_short_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/short/{{sample}}/assembly.fasta"
    output:
        idx = f"{RESULTS_DIR}/indices/single/short/contigs/{{sample}}.1.bt2"
    conda:
        "envs/mapping.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/indices/single/short/contigs
        bowtie2-build {input.contigs} {RESULTS_DIR}/indices/single/short/contigs/{wildcards.sample}
        """

rule index_short_coassembly:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta"
    output:
        idxprefix = f"{RESULTS_DIR}/indices/coassembly/short/contigs/contigs.1.bt2"
    conda:
        "envs/mapping.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/indices/coassembly/short/contigs
        bowtie2-build {input.contigs} \
          {RESULTS_DIR}/indices/coassembly/short/contigs/contigs
        """

rule index_hybrid_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/hybrid/{{sample}}/assembly.fasta"
    output:
        idx = f"{RESULTS_DIR}/indices/single/hybrid/contigs/{{sample}}.1.bt2"
    conda:
        "envs/mapping.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/indices/single/hybrid/contigs
        bowtie2-build {input.contigs} {RESULTS_DIR}/indices/single/hybrid/contigs/{wildcards.sample}
        """

############################################
# 5.2 Mapping
############################################

# ---- SHORT READS ----

rule map_short_single:
    input:
        idx = f"{RESULTS_DIR}/indices/single/short/contigs/{{sample}}.1.bt2",
        r1 = SHORT_FINAL_R1,
        r2 = SHORT_FINAL_R2
    output:
        bam = f"{RESULTS_DIR}/mapping/single/short/{{sample}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        """
        PREFIX=$(basename {input.idx} .1.bt2)
        bowtie2 -x {RESULTS_DIR}/indices/single/short/contigs/$PREFIX -1 {input.r1} -2 {input.r2} -p {threads} |
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule map_short_coassembly:
    input:
        idxprefix = f"{RESULTS_DIR}/indices/coassembly/short/contigs/contigs.1.bt2",
        r1 = SHORT_FINAL_R1,
        r2 = SHORT_FINAL_R2
    output:
        bam = f"{RESULTS_DIR}/mapping/coassembly/short/{{sample}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        """
        PREFIX=$(basename {input.idxprefix} .1.bt2)
        bowtie2 -x results/debug/indices/coassembly/short/contigs/$PREFIX \
          -1 {input.r1} -2 {input.r2} -p {threads} |
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """


rule map_short_multi:
    input:
        idx = f"{RESULTS_DIR}/indices/single/short/contigs/{{sample}}.1.bt2",
        r1 = lambda wc: SHORT_FINAL_R1,
        r2 = lambda wc: SHORT_FINAL_R2
    output:
        bam = f"{RESULTS_DIR}/mapping/multi/short/{{sample}}/{{other}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        """
        PREFIX=$(basename {input.idx} .1.bt2)
        bowtie2 -x {RESULTS_DIR}/indices/single/short/contigs/$PREFIX -1 {input.r1} -2 {input.r2} -p {threads} |
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# ---- LONG READS ----

rule map_long_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/long/{{sample}}/assembly.fasta",
        reads = LONG_FINAL
    output:
        bam = f"{RESULTS_DIR}/mapping/single/long/{{sample}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        """
        minimap2 -ax map-hifi {input.contigs} {input.reads} -t {threads} |
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule map_long_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/long/{{sample}}/assembly.fasta",
        reads = lambda wc: LONG_FINAL.format(sample=wc.other)
    output:
        bam = f"{RESULTS_DIR}/mapping/multi/long/{{sample}}/{{other}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        """
        minimap2 -ax map-hifi {input.contigs} {input.reads} -t {threads} |
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# ---- HYBRID ----

rule map_hybrid_single:
    input:
        idx = f"{RESULTS_DIR}/indices/single/hybrid/contigs/{{sample}}.1.bt2",
        contigs = f"{RESULTS_DIR}/assemblies/single/hybrid/{{sample}}/assembly.fasta",
        r1 = SHORT_FINAL_R1,
        r2 = SHORT_FINAL_R2,
        long = LONG_FINAL
    output:
        bam = f"{RESULTS_DIR}/mapping/single/hybrid/{{sample}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        """
        PREFIX=$(basename {input.idx} .1.bt2)
        bowtie2 -x {RESULTS_DIR}/indices/single/hybrid/contigs/$PREFIX -1 {input.r1} -2 {input.r2} -p {threads} |
        samtools sort -@ {threads} -o short_{wildcards.sample}.bam

        minimap2 -ax map-hifi {input.contigs} {input.long} -t {threads} |
        samtools sort -@ {threads} -o long_{wildcards.sample}.bam

        samtools merge -@ {threads} {output.bam} short_{wildcards.sample}.bam long_{wildcards.sample}.bam
        samtools index {output.bam}

        rm short_{wildcards.sample}.bam long_{wildcards.sample}.bam
        """

rule map_hybrid_multi:
    input:
        idx = f"{RESULTS_DIR}/indices/single/hybrid/contigs/{{sample}}.1.bt2",
        contigs = f"{RESULTS_DIR}/assemblies/single/hybrid/{{sample}}/assembly.fasta",
        r1 = lambda wc: SHORT_FINAL_R1,
        r2 = lambda wc: SHORT_FINAL_R2,
        long = lambda wc: LONG_FINAL
    output:
        bam = f"{RESULTS_DIR}/mapping/multi/hybrid/{{sample}}/{{other}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        """
        PREFIX=$(basename {input.idx} .1.bt2)
        bowtie2 -x {RESULTS_DIR}/indices/single/hybrid/contigs/$PREFIX -1 {input.r1} -2 {input.r2} -p {threads} |
        samtools sort -@ {threads} -o short_{wildcards.sample}_{wildcards.other}.bam

        minimap2 -ax map-hifi {input.contigs} {input.long} -t {threads} |
        samtools sort -@ {threads} -o long_{wildcards.sample}_{wildcards.other}.bam

        samtools merge -@ {threads} {output.bam} \
          short_{wildcards.sample}_{wildcards.other}.bam \
          long_{wildcards.sample}_{wildcards.other}.bam
        samtools index {output.bam}

        rm short_{wildcards.sample}_{wildcards.other}.bam \
           long_{wildcards.sample}_{wildcards.other}.bam
        """

############################################
# 5.3 Depth calculation
############################################

rule depth_coassembly:
    input:
        expand(f"{RESULTS_DIR}/mapping/coassembly/short/{{sample}}.sorted.bam", sample=SAMPLES)
    output:
        depth = f"{RESULTS_DIR}/depth/coassembly/depth.txt"
    conda:
        "envs/binning.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input}
        """

rule depth_single:
    input:
        bam = f"{RESULTS_DIR}/mapping/single/{{assembly_type}}/{{sample}}.sorted.bam"
    output:
        depth = f"{RESULTS_DIR}/depth/single/{{assembly_type}}/{{sample}}/depth.txt"
    conda:
        "envs/binning.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        """

rule depth_multi:
    input:
        bams = lambda wc: expand(
            f"{RESULTS_DIR}/mapping/multi/{wc.assembly_type}/{wc.sample}/{{other}}.sorted.bam",
            other=SAMPLES
        )
    output:
        depth = f"{RESULTS_DIR}/depth/multi/{{assembly_type}}/{{sample}}/depth.txt"
    conda:
        "envs/binning.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bams}
        """

############################################
# 5.4 MetaBAT2
############################################

rule metabat2_coassembly:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta",
        depth = f"{RESULTS_DIR}/depth/coassembly/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/metabat2/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/bins/coassembly/metabat2
        metabat2 -t {threads} -i {input.contigs} -a {input.depth} \
          -o {RESULTS_DIR}/bins/coassembly/metabat2/bin
        touch {output}
        """

rule metabat2_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        depth = f"{RESULTS_DIR}/depth/single/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/single/metabat2/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/bins/single/metabat2/{wildcards.assembly_type}/{wildcards.sample}
        metabat2 -t {threads} -i {input.contigs} -a {input.depth} \
          -o {RESULTS_DIR}/bins/single/metabat2/{wildcards.assembly_type}/{wildcards.sample}/bin
        touch {output}
        """

rule metabat2_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        depth = f"{RESULTS_DIR}/depth/multi/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/multi/metabat2/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/bins/multi/metabat2/{wildcards.assembly_type}/{wildcards.sample}
        metabat2 -t {threads} \
          -i {input.contigs} \
          -a {input.depth} \
          -o {RESULTS_DIR}/bins/multi/metabat2/{wildcards.assembly_type}/{wildcards.sample}/bin
        touch {output}
        """

############################################
# 5.5 MaxBin 2
############################################

rule maxbin2_coassembly:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta",
        depth   = f"{RESULTS_DIR}/depth/coassembly/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/maxbin2/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/bins/coassembly/maxbin2
        run_MaxBin.pl \
          -contig {input.contigs} \
          -abund {input.depth} \
          -out {RESULTS_DIR}/bins/coassembly/maxbin2/bin \
          -thread 2
        touch {output}
        """

rule maxbin2_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        depth   = f"{RESULTS_DIR}/depth/single/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/single/maxbin2/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    log:
        f"logs/{DATASET}/binning/maxbin2/single/{{assembly_type}}/{{sample}}.log"
    shell:
        """
        mkdir -p $(dirname {log})
        mkdir -p {RESULTS_DIR}/bins/single/maxbin2/{wildcards.assembly_type}/{wildcards.sample}

        run_MaxBin.pl \
          -contig {input.contigs} \
          -abund {input.depth} \
          -out {RESULTS_DIR}/bins/single/maxbin2/{wildcards.assembly_type}/{wildcards.sample}/bin \
          -thread 2 \
          > {log} 2>&1

        touch {output}
        """

rule maxbin2_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        depth   = f"{RESULTS_DIR}/depth/multi/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/multi/maxbin2/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        mkdir -p {RESULTS_DIR}/bins/multi/maxbin2/{wildcards.assembly_type}/{wildcards.sample}
        echo {input.depth}
        run_MaxBin.pl \
          -contig {input.contigs} \
          -abund {input.depth} \
          -out {RESULTS_DIR}/bins/multi/maxbin2/{wildcards.assembly_type}/{wildcards.sample}/bin \
          -thread 2
        touch {output}
        """

############################################
# 5.6 CONCOCT
############################################

rule concoct_coassembly:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta",
        bams = expand(
            f"{RESULTS_DIR}/mapping/coassembly/short/{{sample}}.sorted.bam",
            sample=SAMPLES
        )
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/concoct/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/coassembly/concoct
        mkdir -p $outdir

        cut_up_fasta.py {input.contigs} \
          -c 10000 -o 0 --merge_last \
          -b $outdir/contigs_10K.bed \
          > $outdir/contigs_10K.fa

        concoct_coverage_table.py \
          $outdir/contigs_10K.bed \
          {input.bams} \
          > $outdir/coverage_table.tsv

        concoct \
          -t {threads} \
          --composition_file $outdir/contigs_10K.fa \
          --coverage_file $outdir/coverage_table.tsv \
          -b $outdir/concoct_output/

        merge_cutup_clustering.py \
          $outdir/concoct_output/clustering_gt1000.csv \
          > $outdir/concoct_output/clustering_merged.csv

        mkdir -p $outdir/concoct_output/fasta_bins
        extract_fasta_bins.py \
          {input.contigs} \
          $outdir/concoct_output/clustering_merged.csv \
          --output_path $outdir/concoct_output/fasta_bins

        touch {output}
        """

rule concoct_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bam = f"{RESULTS_DIR}/mapping/single/{{assembly_type}}/{{sample}}.sorted.bam"
    output:
        touch(f"{RESULTS_DIR}/bins/single/concoct/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/single/concoct/{wildcards.assembly_type}/{wildcards.sample}
        mkdir -p $outdir

        cut_up_fasta.py {input.contigs} \
          -c 10000 -o 0 --merge_last \
          -b $outdir/contigs_10K.bed \
          > $outdir/contigs_10K.fa

        concoct_coverage_table.py \
          $outdir/contigs_10K.bed \
          {input.bam} \
          > $outdir/coverage_table.tsv

        concoct \
          -t {threads} \
          --composition_file $outdir/contigs_10K.fa \
          --coverage_file $outdir/coverage_table.tsv \
          -b $outdir/concoct_output/

        merge_cutup_clustering.py \
          $outdir/concoct_output/clustering_gt1000.csv \
          > $outdir/concoct_output/clustering_merged.csv

        mkdir -p $outdir/concoct_output/fasta_bins
        extract_fasta_bins.py \
          {input.contigs} \
          $outdir/concoct_output/clustering_merged.csv \
          --output_path $outdir/concoct_output/fasta_bins

        touch {output}
        """

rule concoct_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bams = lambda wc: expand(
            f"{RESULTS_DIR}/mapping/multi/{wc.assembly_type}/{wc.sample}/{{other}}.sorted.bam",
            other=SAMPLES
        )
    output:
        touch(f"{RESULTS_DIR}/bins/multi/concoct/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/multi/concoct/{wildcards.assembly_type}/{wildcards.sample}
        mkdir -p $outdir

        cut_up_fasta.py {input.contigs} \
          -c 10000 -o 0 --merge_last \
          -b $outdir/contigs_10K.bed \
          > $outdir/contigs_10K.fa

        concoct_coverage_table.py \
          $outdir/contigs_10K.bed \
          {input.bams} \
          > $outdir/coverage_table.tsv

        concoct \
          -t {threads} \
          --composition_file $outdir/contigs_10K.fa \
          --coverage_file $outdir/coverage_table.tsv \
          -b $outdir/concoct_output/

        merge_cutup_clustering.py \
          $outdir/concoct_output/clustering_gt1000.csv \
          > $outdir/concoct_output/clustering_merged.csv

        mkdir -p $outdir/concoct_output/fasta_bins
        extract_fasta_bins.py \
          {input.contigs} \
          $outdir/concoct_output/clustering_merged.csv \
          --output_path $outdir/concoct_output/fasta_bins

        touch {output}
        """

############################################
# 5.7. VAMB
############################################

rule vamb_coassembly:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta",
        depth   = f"{RESULTS_DIR}/depth/coassembly/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/vamb/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/coassembly/vamb
        mkdir -p $outdir

        vamb \
          --outdir $outdir \
          --fasta {input.contigs} \
          --jgi {input.depth} \
          --minfasta 200000 \
          -m 2000 \
          --threads {threads}

        touch {output}
        """

rule vamb_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        depth   = f"{RESULTS_DIR}/depth/single/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/single/vamb/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/single/vamb/{wildcards.assembly_type}/{wildcards.sample}
        mkdir -p $outdir

        vamb \
          --outdir $outdir \
          --fasta {input.contigs} \
          --jgi {input.depth} \
          --minfasta 200000 \
          -m 2000 \
          --threads {threads}

        touch {output}
        """

rule vamb_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        depth   = f"{RESULTS_DIR}/depth/multi/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/multi/vamb/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/multi/vamb/{wildcards.assembly_type}/{wildcards.sample}
        mkdir -p $outdir

        vamb \
          --outdir $outdir \
          --fasta {input.contigs} \
          --jgi {input.depth} \
          --minfasta 200000 \
          -m 2000 \
          --threads {threads}

        touch {output}
        """

############################################
# 5.8. MetaDecoder
############################################

rule metadecoder_coassembly:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta",
        sams = expand(
            f"{RESULTS_DIR}/mapping/coassembly/short/{{sample}}.sorted.bam",
            sample=SAMPLES
        )

    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/metadecoder/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_metadecoder.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/coassembly/metadecoder
        mkdir -p $outdir
        cd $outdir

        metadecoder coverage \
          --threads {threads} \
          -s {input.sams} \
          -o METADECODER_gsa.COVERAGE

        metadecoder seed \
          --threads {threads} \
          -f {input.contigs} \
          -o METADECODER_gsa.SEED

        metadecoder cluster \
          -f {input.contigs} \
          -c METADECODER_gsa.COVERAGE \
          -s METADECODER_gsa.SEED \
          -o METADECODER_coassembly

        touch {output}
        """

rule metadecoder_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        sam = f"{RESULTS_DIR}/mapping/single/{{assembly_type}}/{{sample}}/{{sample}}.sorted.bam"
    output:
        touch(f"{RESULTS_DIR}/bins/single/metadecoder/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_metadecoder.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/single/metadecoder/{wildcards.assembly_type}/{wildcards.sample}
        mkdir -p $outdir
        cd $outdir

        metadecoder coverage \
          --threads {threads} \
          -s {input.sam} \
          -o METADECODER_gsa.COVERAGE

        metadecoder seed \
          --threads {threads} \
          -f {input.contigs} \
          -o METADECODER_gsa.SEED

        metadecoder cluster \
          -f {input.contigs} \
          -c METADECODER_gsa.COVERAGE \
          -s METADECODER_gsa.SEED \
          -o METADECODER_{wildcards.sample}

        touch {output}
        """

rule metadecoder_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        sams = lambda wc: expand(
            f"{RESULTS_DIR}/mapping/multi/{wc.assembly_type}/{wc.sample}/{{rep}}.sorted.bam",
            rep=SAMPLES
        )

    output:
        touch(f"{RESULTS_DIR}/bins/multi/metadecoder/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_metadecoder.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/multi/metadecoder/{wildcards.assembly_type}/{wildcards.sample}
        mkdir -p $outdir
        cd $outdir

        metadecoder coverage \
          --threads {threads} \
          -s {input.sams} \
          -o METADECODER_gsa.COVERAGE

        metadecoder seed \
          --threads {threads} \
          -f {input.contigs} \
          -o METADECODER_gsa.SEED

        metadecoder cluster \
          -f {input.contigs} \
          -c METADECODER_gsa.COVERAGE \
          -s METADECODER_gsa.SEED \
          -o METADECODER_{wildcards.sample}

        touch {output}
        """

############################################
# 5.9. BINNY
############################################

rule binny_coassembly:
    input:
        config = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta"
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/binny/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/coassembly/binny
        mkdir -p $outdir

        binny \
          -l \
          -n marine_result \
          -r \
          -t {threads} \
          {input.config}

        touch {output}
        """

rule binny_single:
    input:
        config = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta"
    output:
        touch(f"{RESULTS_DIR}/bins/single/binny/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/single/binny/{wildcards.sample}
        mkdir -p $outdir

        binny \
          -l \
          -n {wildcards.sample}_single_result \
          -r \
          -t {threads} \
          {input.config}

        touch {output}
        """

rule binny_multi:
    input:
        config = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta"
    output:
        touch(f"{RESULTS_DIR}/bins/multi/binny/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/multi/binny/{wildcards.sample}
        mkdir -p $outdir

        binny \
          -l \
          -n {wildcards.sample}_multi_result \
          -r \
          -t {threads} \
          {input.config}

        touch {output}
        """

############################################
# 5.10. MetaBinner
############################################

# TODO metabinner_path in shell command should probably be: $CONDA_PREFIX/bin/MetaBinner

rule metabinner_coassembly:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta",
        coverage = f"{RESULTS_DIR}/depth/coassembly/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/metabinner/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_metabinner.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/coassembly/metabinner
        mkdir -p $outdir

        metabinner_path=xx/MetaBinner
        python xx/MetaBinner/scripts/gen_kmer.py {input.contigs} 1000 4 $outdir/marine_kmer.tsv

        contig_file={input.contigs}
        output_dir=$outdir/output
        coverage_profiles={input.coverage}
        kmer_profile=$outdir/marine_kmer.tsv

        bash xx/MetaBinner/scripts/run_metabinner.sh \
          -t {threads} \
          -a $contig_file \
          -o $output_dir \
          -d $coverage_profiles \
          -k $kmer_profile \
          -p $metabinner_path

        touch {output}
        """

rule metabinner_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        coverage = f"{RESULTS_DIR}/depth/single/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/single/metabinner/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_metabinner.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/single/metabinner/{{assembly_type}}/{{wildcards.sample}}
        mkdir -p $outdir

        metabinner_path=xx/MetaBinner
        python xx/MetaBinner/scripts/gen_kmer.py {input.contigs} 1000 4 $outdir/{{wildcards.sample}}_kmer.tsv

        contig_file={input.contigs}
        output_dir=$outdir/output
        coverage_profiles={input.coverage}
        kmer_profile=$outdir/{{wildcards.sample}}_kmer.tsv

        bash xx/MetaBinner/scripts/run_metabinner.sh \
          -t {threads} \
          -a $contig_file \
          -o $output_dir \
          -d $coverage_profiles \
          -k $kmer_profile \
          -p $metabinner_path

        touch {output}
        """

rule metabinner_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        coverage = f"{RESULTS_DIR}/depth/multi/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/multi/metabinner/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_metabinner.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/multi/metabinner/{{assembly_type}}/{{wildcards.sample}}
        mkdir -p $outdir

        metabinner_path=xx/MetaBinner
        python xx/MetaBinner/scripts/gen_kmer.py {input.contigs} 1000 4 $outdir/{{wildcards.sample}}_kmer.tsv

        contig_file={input.contigs}
        output_dir=$outdir/output
        coverage_profiles={input.coverage}
        kmer_profile=$outdir/{{wildcards.sample}}_kmer.tsv

        bash xx/MetaBinner/scripts/run_metabinner.sh \
          -t {threads} \
          -a $contig_file \
          -o $output_dir \
          -d $coverage_profiles \
          -k $kmer_profile \
          -p $metabinner_path

        touch {output}
        """

############################################
# 5.11. SemiBin2
############################################

rule semibin2_coassembly:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta",
        bams = expand(f"{RESULTS_DIR}/mapping/coassembly/short/{{sample}}.sorted.bam", sample=SAMPLES)
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/semibin2/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    params:
        bam_list=lambda wc, input: ' '.join(input.bams)
    shell:
        """
        outdir={RESULTS_DIR}/bins/coassembly/semibin2
        mkdir -p $outdir

        SemiBin2 single_easy_bin \
          -t {threads} \
          -i {input.contigs} \
          -b bam_list \
          -o $outdir/output \
          --compression none

        touch {output}
        """

rule semibin2_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bam = f"{RESULTS_DIR}/mapping/single/{{assembly_type}}/{{sample}}/{{sample}}.sorted.bam"
    output:
        touch(f"{RESULTS_DIR}/bins/single/semibin2/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/single/semibin2/{{assembly_type}}/{{wildcards.sample}}
        mkdir -p $outdir

        SemiBin2 single_easy_bin \
          -t {threads} \
          -i {input.contigs} \
          -b {input.bam} \
          -o $outdir/output \
          --compression none

        touch {output}
        """

rule semibin2_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bams = lambda wc: expand(
            f"{RESULTS_DIR}/mapping/multi/{wc.assembly_type}/{wc.sample}/{{other}}.sorted.bam",
            other=SAMPLES
        )    
    output:
        touch(f"{RESULTS_DIR}/bins/multi/semibin2/{{assembly_type}}/{{sample}}/bins.done")
    threads: 
        config["threads"]
    conda:
        "envs/binning.yaml"
    params:
        bam_list=lambda wc, input: ' '.join(input.bams)
    shell:
        """
        outdir={RESULTS_DIR}/bins/multi/semibin2/{{assembly_type}}/{{wildcards.sample}}
        mkdir -p $outdir

        SemiBin2 single_easy_bin \
          -t {threads} \
          -i {input.contigs} \
          -b bam_list \
          -o $outdir/output \
          --compression none

        touch {output}
        """

############################################
# 5.12. COMEBin
############################################

rule comebin_coassembly:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta",
        bams_dir = expand(f"{RESULTS_DIR}/mapping/coassembly/short/{{sample}}.sorted.bam", sample=SAMPLES)
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/comebin/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_comebin.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/coassembly/comebin
        mkdir -p $outdir

        run_comebin.sh \
          -t {threads} \
          -a {input.contigs} \
          -o $outdir \
          -p {input.bams_dir}

        touch {output}
        """

rule comebin_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bam = f"{RESULTS_DIR}/mapping/single/{{assembly_type}}/{{sample}}/{{sample}}.sorted.bam"
    output:
        touch(f"{RESULTS_DIR}/bins/single/comebin/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_comebin.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/single/comebin/{{assembly_type}}/{{wildcards.sample}}
        mkdir -p $outdir

        run_comebin.sh \
          -t {threads} \
          -a {input.contigs} \
          -o $outdir \
          -p {input.bam}

        touch {output}
        """

rule comebin_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bams = lambda wc: expand(
            f"{RESULTS_DIR}/mapping/multi/{wc.assembly_type}/{wc.sample}/{{other}}.sorted.bam",
            other=SAMPLES
        ) 
    output:
        touch(f"{RESULTS_DIR}/bins/multi/comebin/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_comebin.yaml"
    shell:
        """
        outdir={RESULTS_DIR}/bins/multi/comebin/{{assembly_type}}/{{wildcards.sample}}
        mkdir -p $outdir

        run_comebin.sh \
          -t {threads} \
          -a {input.contigs} \
          -o $outdir \
          -p {input.bams}

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
