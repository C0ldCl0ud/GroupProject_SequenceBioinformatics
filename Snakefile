wildcard_constraints: 
    sample = "S\\d+"

import os
from glob import glob

configfile: "config/config.yaml"

# dataset is passed via --config dataset=dataset1
DATASET = config["dataset"]
DATASET_DIR = f"data/{DATASET}"
RESULTS_DIR = f"results/{DATASET}"

USE_SCRATCH = config.get("use_scratch", False)

if USE_SCRATCH:
    SCRATCH_BASE = config.get("scratch_base", "/beegfs/HPCscratch")
    USER = os.environ.get("USER", "unknown")
    PROJECT = config.get("dataset", "GroupProject")
    SCRATCH = f"{SCRATCH_BASE}/{USER}/{PROJECT}"
    SCRATCH_MAP = f"{SCRATCH}/mapping"
else:
    SCRATCH = ""
    SCRATCH_MAP = f"{RESULTS_DIR}/mapping"

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

DEBUG = config.get("debug", False)

VAMB_MINFASTA = (
    config["vamb"]["debug"]["minfasta"]
    if DEBUG else
    config["vamb"]["prod"]["minfasta"]
)

VAMB_MINCONTIGS = (
    config["vamb"]["debug"]["mincontigs"]
    if DEBUG else
    config["vamb"]["prod"]["mincontigs"]
)

VAMB_BATCHSIZE = (
    config["vamb"]["debug"]["batchsize"]
    if DEBUG else
    config["vamb"]["prod"]["batchsize"]
)

############################################
# Rule order / final targets
############################################

rule all:
    input:
        [
        #expand(f"{RESULTS_DIR}/raw/short/{{sample}}", sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/raw/long/{{sample}}", sample=SAMPLES)

        #expand(f"{RESULTS_DIR}/fastq/short/{{sample}}_R1.fq.gz", sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/fastq/short/{{sample}}_R2.fq.gz", sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/fastq/long/{{sample}}.fq.gz", sample=SAMPLES),

        #expand(f"{RESULTS_DIR}/preprocess/short/{{sample}}_R1.fastp.fq.gz", sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/preprocess/short/{{sample}}_R2.fastp.fq.gz", sample=SAMPLES)
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
        #expand(f"{RESULTS_DIR}/bins/coassembly/{{tool}}/bins.done",tool=config["binning_tools"]),
        #expand(f"{RESULTS_DIR}/bins/single/{{tool}}/{{assembly_type}}/{{sample}}/bins.done", tool=config["binning_tools"], assembly_type=["short", "long", "hybrid"], sample=SAMPLES),
        #expand(f"{RESULTS_DIR}/bins/multi/{{tool}}/{{assembly_type}}/{{sample}}/bins.done", tool=config["binning_tools"], sample=SAMPLES, assembly_type=["short", "long", "hybrid"])

        # evaluation
        # Single-sample eval
        expand(f"{RESULTS_DIR}/eval/single/{{tool}}/{{assembly_type}}/{{sample}}/eval_{{eval_type}}.done",
                tool=config["binning_tools"],
                assembly_type=["short","long","hybrid"],
                sample=SAMPLES,
                eval_type=["comp_cont", "tRNA", "rRNA"]),
                #eval_type=["comp_cont"]),
                #eval_type=["rRNA"]),
        # Multi-sample eval
        expand(f"{RESULTS_DIR}/eval/multi/{{tool}}/{{assembly_type}}/{{sample}}/eval_{{eval_type}}.done",
                tool=config["binning_tools"],
                assembly_type=["short","long","hybrid"],
                sample=SAMPLES,
                eval_type=["comp_cont", "tRNA", "rRNA"]),
                #eval_type=["comp_cont"]),
                #eval_type=["rRNA"]),
        # Coassembly eval
        expand(f"{RESULTS_DIR}/eval/coassembly/{{tool}}/eval_{{eval_type}}.done",
                tool=config["binning_tools"],
                eval_type=["comp_cont", "tRNA", "rRNA"]),
                #eval_type=["comp_cont"]),
                #eval_type=["rRNA"])
        ]

############################################
# 1. Download SRA / dump files
############################################

rule download_short_reads:
    output:
        sra = directory(f"{RESULTS_DIR}/raw/short/{{sample}}")
    params:
        acc = lambda wc: SHORT_ACC[wc.sample],
        scratch = SCRATCH
    conda:
        "envs/download.yaml"
    log:
        f"logs/{DATASET}/download/short/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail

        if [ -n "{params.scratch}" ]; then
            mkdir -p "{params.scratch}"
            chmod 700 "{params.scratch}"
            TMPDIR=$(mktemp -d "{params.scratch}/prefetch_{wildcards.sample}_XXXX")
        else
            TMPDIR=$(mktemp -d)
        fi

        trap "rm -rf $TMPDIR" EXIT

        prefetch {params.acc} -X 40G -O "$TMPDIR" > {log} 2>&1

        mkdir -p $(dirname {output.sra})
        mv "$TMPDIR" "{output.sra}"
        """

rule download_long_reads:
    output:
        sra = directory(f"{RESULTS_DIR}/raw/long/{{sample}}")
    params:
        acc = lambda wc: LONG_ACC[wc.sample],
        scratch = SCRATCH
    conda:
        "envs/download.yaml"
    log:
        f"logs/{DATASET}/download/long/{{sample}}.log"
    resources:
        disk_io=1
    shell:
        r"""
        set -euo pipefail

        if [ -n "{params.scratch}" ]; then
            mkdir -p "{params.scratch}"
            chmod 700 "{params.scratch}"
            TMPDIR=$(mktemp -d "{params.scratch}/prefetch_{wildcards.sample}_XXXX")
        else
            TMPDIR=$(mktemp -d)
        fi

        trap "rm -rf $TMPDIR" EXIT

        prefetch {params.acc} -X 40G -O "$TMPDIR" > {log} 2>&1

        mkdir -p $(dirname {output.sra})
        mv "$TMPDIR" "{output.sra}"
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
        seed = config.get("debug_seed", 42),
        scratch = SCRATCH
    conda:
        "envs/download.yaml"
    log:
        f"logs/{DATASET}/fastq/short/{{sample}}.log"
    resources:
        disk_io=1
    shell:
        r"""
        set -euo pipefail

        SRA="{input.sra_dir}/{params.acc}/{params.acc}.sra"

        # --- scratch handling ---
        if [ -n "{params.scratch}" ] && [ -d "{params.scratch}" ] && [ -w "{params.scratch}" ]; then
            mkdir -p "{params.scratch}"
            chmod 700 "{params.scratch}"
            TMPDIR=$(mktemp -d "{params.scratch}/fastq_short_{wildcards.sample}_XXXX")
        else
            TMPDIR=$(mktemp -d)
        fi

        trap "rm -rf $TMPDIR" EXIT

        R1_UNZ="$TMPDIR/{params.acc}_1.fastq"
        R2_UNZ="$TMPDIR/{params.acc}_2.fastq"

        # --- extract FASTQ into scratch ---
        fasterq-dump --split-files "$SRA" -O "$TMPDIR" >> {log} 2>&1

        mkdir -p $(dirname "{output.r1}")

        if [ "{params.debug}" = "True" ]; then
            seqtk sample -s{params.seed} "$R1_UNZ" {params.n} | gzip > "{output.r1}"
            seqtk sample -s{params.seed} "$R2_UNZ" {params.n} | gzip > "{output.r2}"
        else
            gzip -c "$R1_UNZ" > "{output.r1}"
            gzip -c "$R2_UNZ" > "{output.r2}"
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
        seed = config.get("debug_seed", 42),
        scratch = SCRATCH
    conda:
        "envs/download.yaml"
    log:
        f"logs/{DATASET}/fastq/long/{{sample}}.log"
    resources:
        disk_io=1
    shell:
        r"""
        set -euo pipefail

        SRA="{input.sra_dir}/{params.acc}/{params.acc}.sra"

        # --- scratch handling ---
        if [ -n "{params.scratch}" ] && [ -d "{params.scratch}" ] && [ -w "{params.scratch}" ]; then
            mkdir -p "{params.scratch}"
            chmod 700 "{params.scratch}"
            TMPDIR=$(mktemp -d "{params.scratch}/fastq_long_{wildcards.sample}_XXXX")
        else
            TMPDIR=$(mktemp -d)
        fi

        trap "rm -rf $TMPDIR" EXIT

        UNZ="$TMPDIR/{params.acc}.fastq"

        # --- extract FASTQ into scratch ---
        fasterq-dump "$SRA" -O "$TMPDIR" >> {log} 2>&1

        mkdir -p $(dirname "{output.fq}")

        if [ "{params.debug}" = "True" ]; then
            seqtk sample -s{params.seed} "$UNZ" {params.n} | gzip > "{output.fq}"
        else
            gzip -c "$UNZ" > "{output.fq}"
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
        r"""
        set -euo pipefail

        TMPDIR="{RESULTS_DIR}/tmp/{wildcards.sample}"
        mkdir -p "$TMPDIR"

        gunzip -c {input.fq} > "$TMPDIR/{wildcards.sample}.fq"

        qcat --trim --detect-middle \
             -f "$TMPDIR/{wildcards.sample}.fq" \
             -o "$TMPDIR/{wildcards.sample}_qcat.fq" \
            > {log} 2>&1

        gzip -c "$TMPDIR/{wildcards.sample}_qcat.fq" > {output.fq}

        rm -rf "$TMPDIR"
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
          --un-conc-gz {RESULTS_DIR}/preprocess/short/{wildcards.sample} \
          -S /dev/null \
          > {log} 2>&1

        mv {RESULTS_DIR}/preprocess/short/{wildcards.sample}.1 {output.r1}
        mv {RESULTS_DIR}/preprocess/short/{wildcards.sample}.2 {output.r2}
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
          --un-gz {RESULTS_DIR}/preprocess/long/{wildcards.sample} \
          -S /dev/null \
          > {log} 2>&1

        mv {RESULTS_DIR}/preprocess/long/{wildcards.sample} {output}
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
    threads: 16
    log:
        f"logs/{DATASET}/assembly/single/long/{{sample}}.flye.log"
    resources:
        ram=10
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
        "containers/operams.simg" # remote server
        #"containers/operams.sif"   # local machine
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
        idxdir = directory(f"{RESULTS_DIR}/indices/coassembly/short/contigs")
    conda:
        "envs/mapping.yaml"
    shell:
        """
        mkdir -p {output.idxdir}
        bowtie2-build {input.contigs} {output.idxdir}/contigs
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
        bam = f"{SCRATCH_MAP}/single/short/{{sample}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.bam})
        rm -f {output.bam} {output.bam}.bai

        # Create UNIQUE temp dir per job
        if [ -n "{SCRATCH}" ]; then
            TMPDIR=$(mktemp -d {SCRATCH}/samtools.{wildcards.sample}.XXXXXX)
        else
            TMPDIR=$(mktemp -d ./tmp_samtools.XXXXXX)
        fi

        # Always clean up
        trap 'rm -rf "$TMPDIR"' EXIT

        PREFIX=$(basename {input.idx} .1.bt2)

        bowtie2 \
          -x {RESULTS_DIR}/indices/single/short/contigs/$PREFIX \
          -1 {input.r1} -2 {input.r2} \
          -p {threads} |
        samtools sort \
          -@ {threads} \
          -T "$TMPDIR/sort" \
          -o {output.bam}

        samtools index {output.bam}
        """


rule map_short_coassembly:
    input:
        idxdir = f"{RESULTS_DIR}/indices/coassembly/short/contigs",
        r1 = SHORT_FINAL_R1,
        r2 = SHORT_FINAL_R2
    output:
        bam = f"{SCRATCH_MAP}/coassembly/short/{{sample}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        """
        bowtie2 -x {input.idxdir}/contigs \
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
        bam = f"{SCRATCH_MAP}/multi/short/{{sample}}/{{other}}.sorted.bam",
        bai = f"{SCRATCH_MAP}/multi/short/{{sample}}/{{other}}.sorted.bam.bai"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.bam})
        rm -f {output.bam} {output.bai}

        # Create UNIQUE temp directory per job
        if [ -n "{SCRATCH}" ]; then
            TMPDIR=$(mktemp -d {SCRATCH}/samtools.{wildcards.sample}.{wildcards.other}.XXXXXX)
        else
            TMPDIR=$(mktemp -d ./tmp_samtools.XXXXXX)
        fi

        # Always clean up temp dir
        trap 'rm -rf "$TMPDIR"' EXIT

        PREFIX=$(basename {input.idx} .1.bt2)

        bowtie2 \
            -x {RESULTS_DIR}/indices/single/short/contigs/$PREFIX \
            -1 {input.r1} \
            -2 {input.r2} \
            -p {threads} |
        samtools sort \
            -@ {threads} \
            -T "$TMPDIR/sort" \
            -o {output.bam}

        samtools index {output.bam}
        """


# ---- LONG READS ----

rule map_long_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/long/{{sample}}/assembly.fasta",
        reads = LONG_FINAL
    output:
        bam = f"{SCRATCH_MAP}/single/long/{{sample}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        """
        # Ensure output directory exists
        mkdir -p $(dirname {output.bam})

        # Remove old partial BAM/index if present
        rm -rf {output.bam} {output.bam}.tmp*

        # Create unique temp dir per job
        if [ -n "{SCRATCH}" ]; then
            TMPDIR=$(mktemp -d {SCRATCH}/samtools.{wildcards.sample}.XXXXXX)
        else
            TMPDIR=$(mktemp -d ./tmp_samtools.{wildcards.sample}.XXXXXX)
        fi

        # Always clean up
        trap 'rm -rf "$TMPDIR"' EXIT

        minimap2 -ax map-hifi {input.contigs} {input.reads} -t {threads} |
        samtools sort -@ {threads} -T "$TMPDIR/{wildcards.sample}.tmp" -o {output.bam}
        samtools index {output.bam}

        samtools index {output.bam}
        """

rule map_long_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/long/{{sample}}/assembly.fasta",
        reads = lambda wc: LONG_FINAL.format(sample=wc.other)
    output:
        bam = f"{SCRATCH_MAP}/multi/long/{{sample}}/{{other}}.sorted.bam",
        bai = f"{SCRATCH_MAP}/multi/long/{{sample}}/{{other}}.sorted.bam.bai"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        r"""
        set -euo pipefail

        mkdir -p $(dirname {output.bam})
        rm -f {output.bam} {output.bai}

        # Create a UNIQUE temp dir per job
        if [ -n "{SCRATCH}" ]; then
            TMPDIR=$(mktemp -d {SCRATCH}/samtools.{wildcards.sample}.{wildcards.other}.XXXXXX)
        else
            TMPDIR=$(mktemp -d ./tmp_samtools.XXXXXX)
        fi

        # Always clean up
        trap 'rm -rf "$TMPDIR"' EXIT

        minimap2 -ax map-hifi \
            {input.contigs} {input.reads} -t {threads} |
        samtools sort -@ {threads} \
            -T "$TMPDIR/sort" \
            -o {output.bam}

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
        bam = f"{SCRATCH_MAP}/single/hybrid/{{sample}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        rm -f {output.bam}

        # Unique temp dir
        if [ -n "{SCRATCH}" ]; then
            TMPDIR=$(mktemp -d {SCRATCH}/samtools.{wildcards.sample}.XXXXXX)
        else
            TMPDIR=$(mktemp -d ./tmp_samtools.{wildcards.sample}.XXXXXX)
        fi
        trap 'rm -rf "$TMPDIR"' EXIT

        PREFIX=$(basename {input.idx} .1.bt2)

        # Short reads
        bowtie2 -x {RESULTS_DIR}/indices/single/hybrid/contigs/$PREFIX \
          -1 {input.r1} -2 {input.r2} -p {threads} |
        samtools sort -@ {threads} -T "$TMPDIR/short" -o "$TMPDIR/short.bam"

        # Long reads
        minimap2 -ax map-hifi {input.contigs} {input.long} -t {threads} |
        samtools sort -@ {threads} -T "$TMPDIR/long" -o "$TMPDIR/long.bam"

        # Merge
        samtools merge -@ {threads} {output.bam} "$TMPDIR/short.bam" "$TMPDIR/long.bam"
        samtools index {output.bam}
        """


rule map_hybrid_multi:
    input:
        idx = f"{RESULTS_DIR}/indices/single/hybrid/contigs/{{sample}}.1.bt2",
        contigs = f"{RESULTS_DIR}/assemblies/single/hybrid/{{sample}}/assembly.fasta",
        r1 = lambda wc: SHORT_FINAL_R1,
        r2 = lambda wc: SHORT_FINAL_R2,
        long = lambda wc: LONG_FINAL
    output:
        bam = f"{SCRATCH_MAP}/multi/hybrid/{{sample}}/{{other}}.sorted.bam"
    threads: config["threads"]
    conda:
        "envs/mapping.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        rm -f {output.bam}

        # Unique temp dir
        if [ -n "{SCRATCH}" ]; then
            TMPDIR=$(mktemp -d {SCRATCH}/samtools.{wildcards.sample}.{wildcards.other}.XXXXXX)
        else
            TMPDIR=$(mktemp -d ./tmp_samtools.{wildcards.sample}.{wildcards.other}.XXXXXX)
        fi
        trap 'rm -rf "$TMPDIR"' EXIT

        PREFIX=$(basename {input.idx} .1.bt2)

        # Short reads
        bowtie2 -x {RESULTS_DIR}/indices/single/hybrid/contigs/$PREFIX \
          -1 {input.r1} -2 {input.r2} -p {threads} |
        samtools sort -@ {threads} -T "$TMPDIR/short" -o "$TMPDIR/short.bam"

        # Long reads
        minimap2 -ax map-hifi {input.contigs} {input.long} -t {threads} |
        samtools sort -@ {threads} -T "$TMPDIR/long" -o "$TMPDIR/long.bam"

        # Merge
        samtools merge -@ {threads} {output.bam} "$TMPDIR/short.bam" "$TMPDIR/long.bam"
        samtools index {output.bam}
        """

############################################
# 5.3 Depth calculation
############################################

rule depth_coassembly:
    input:
        expand(f"{SCRATCH_MAP}/coassembly/short/{{sample}}.sorted.bam", sample=SAMPLES)
    output:
        depth = f"{RESULTS_DIR}/depth/coassembly/depth.txt"
    conda:
        "envs/binning.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input}
        """
rule maxbin_abundance_coassembly:
    input:
        depth = f"{RESULTS_DIR}/depth/coassembly/depth.txt"
    output:
        maxbin = f"{RESULTS_DIR}/depth/coassembly/depth.maxbin2.txt"
    shell:
        """
        awk '
            NR==1 {{
            if ($3 ~ /^[0-9.]+$/) print $1 "\t" $3;
            next
            }}
        {{ print $1 "\t" $3 }}
        ' {input.depth} > {output.maxbin}
        """

rule depth_single:
    input:
        bam = f"{SCRATCH_MAP}/single/{{assembly_type}}/{{sample}}.sorted.bam"
    output:
        depth = f"{RESULTS_DIR}/depth/single/{{assembly_type}}/{{sample}}/depth.txt"
    conda:
        "envs/binning.yaml"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        """

rule maxbin_abundance_single:
    input:
        depth = f"{RESULTS_DIR}/depth/single/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        maxbin = f"{RESULTS_DIR}/depth/single/{{assembly_type}}/{{sample}}/depth.maxbin2.txt"
    shell:
        """
        awk '
        NR==1 {{
            if ($3 ~ /^[0-9.]+$/) print $1 "\t" $3;
            next
        }}
        {{ print $1 "\t" $3 }}
        ' {input.depth} > {output.maxbin}
        """

rule depth_multi:
    input:
        bams = lambda wc: expand(
            f"{SCRATCH_MAP}/multi/{wc.assembly_type}/{wc.sample}/{{other}}.sorted.bam",
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

rule maxbin_abundance_multie:
    input:
        depth = f"{RESULTS_DIR}/depth/multi/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        maxbin = f"{RESULTS_DIR}/depth/multi/{{assembly_type}}/{{sample}}/depth.maxbin2.txt"
    shell:
        """
        awk '
        NR==1 {{
            if ($3 ~ /^[0-9.]+$/) print $1 "\t" $3;
            next
        }}
        {{ print $1 "\t" $3 }}
        ' {input.depth} > {output.maxbin}
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
        depth   = f"{RESULTS_DIR}/depth/coassembly/depth.maxbin2.txt"
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
        depth   = f"{RESULTS_DIR}/depth/single/{{assembly_type}}/{{sample}}/depth.maxbin2.txt"
    output:
        touch(f"{RESULTS_DIR}/bins/single/maxbin2/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    log:
        f"logs/{DATASET}/bins/maxbin2/single/{{assembly_type}}/{{sample}}.log"
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
        depth   = f"{RESULTS_DIR}/depth/multi/{{assembly_type}}/{{sample}}/depth.maxbin2.txt"
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
            f"{SCRATCH_MAP}/coassembly/short/{{sample}}.sorted.bam",
            sample=SAMPLES
        )
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/concoct/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR={RESULTS_DIR}/bins/coassembly/concoct
        SCRATCHDIR={SCRATCH}/concoct/coassembly

        mkdir -p $SCRATCHDIR $OUTDIR

        cut_up_fasta.py {input.contigs} \
          -c 10000 -o 0 --merge_last \
          -b $SCRATCHDIR/contigs_10K.bed \
          > $SCRATCHDIR/contigs_10K.fa

        concoct_coverage_table.py \
          $SCRATCHDIR/contigs_10K.bed \
          {input.bams} \
          > $SCRATCHDIR/coverage_table.tsv

        concoct \
          -t {threads} \
          --composition_file $SCRATCHDIR/contigs_10K.fa \
          --coverage_file $SCRATCHDIR/coverage_table.tsv \
          -b $SCRATCHDIR/concoct_output/

        merge_cutup_clustering.py \
          $SCRATCHDIR/concoct_output/clustering_gt1000.csv \
          > $SCRATCHDIR/concoct_output/clustering_merged.csv

        mkdir -p $SCRATCHDIR/concoct_output/fasta_bins
        extract_fasta_bins.py \
          {input.contigs} \
          $SCRATCHDIR/concoct_output/clustering_merged.csv \
          --output_path $SCRATCHDIR/concoct_output/fasta_bins

        rsync -a $SCRATCHDIR/concoct_output/fasta_bins/ $OUTDIR

        touch {output}
        rm -rf $SCRATCHDIR
        """


rule concoct_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bam = f"{SCRATCH_MAP}/single/{{assembly_type}}/{{sample}}.sorted.bam"
    output:
        touch(f"{RESULTS_DIR}/bins/single/concoct/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR={RESULTS_DIR}/bins/single/concoct/{wildcards.assembly_type}/{wildcards.sample}
        SCRATCHDIR={SCRATCH}/concoct/single/{wildcards.assembly_type}/{wildcards.sample}

        mkdir -p $SCRATCHDIR $OUTDIR

        cut_up_fasta.py {input.contigs} \
          -c 10000 -o 0 --merge_last \
          -b $SCRATCHDIR/contigs_10K.bed \
          > $SCRATCHDIR/contigs_10K.fa

        concoct_coverage_table.py \
          $SCRATCHDIR/contigs_10K.bed \
          {input.bam} \
          > $SCRATCHDIR/coverage_table.tsv

        concoct \
          -t {threads} \
          --composition_file $SCRATCHDIR/contigs_10K.fa \
          --coverage_file $SCRATCHDIR/coverage_table.tsv \
          -b $SCRATCHDIR/concoct_output/

        merge_cutup_clustering.py \
          $SCRATCHDIR/concoct_output/clustering_gt1000.csv \
          > $SCRATCHDIR/concoct_output/clustering_merged.csv

        mkdir -p $SCRATCHDIR/concoct_output/fasta_bins
        extract_fasta_bins.py \
          {input.contigs} \
          $SCRATCHDIR/concoct_output/clustering_merged.csv \
          --output_path $SCRATCHDIR/concoct_output/fasta_bins

        rsync -a $SCRATCHDIR/concoct_output/fasta_bins/ $OUTDIR

        touch {output}
        rm -rf $SCRATCHDIR
        """


rule concoct_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bams = lambda wc: expand(
            f"{SCRATCH_MAP}/multi/{wc.assembly_type}/{wc.sample}/{{other}}.sorted.bam",
            other=SAMPLES
        )
    output:
        touch(f"{RESULTS_DIR}/bins/multi/concoct/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        r"""
        set -euo pipefail

        OUTDIR={RESULTS_DIR}/bins/multi/concoct/{wildcards.assembly_type}/{wildcards.sample}
        SCRATCHDIR={SCRATCH}/concoct/multi/{wildcards.assembly_type}/{wildcards.sample}

        mkdir -p $SCRATCHDIR $OUTDIR

        cut_up_fasta.py {input.contigs} \
          -c 10000 -o 0 --merge_last \
          -b $SCRATCHDIR/contigs_10K.bed \
          > $SCRATCHDIR/contigs_10K.fa

        concoct_coverage_table.py \
          $SCRATCHDIR/contigs_10K.bed \
          {input.bams} \
          > $SCRATCHDIR/coverage_table.tsv

        concoct \
          -t {threads} \
          --composition_file $SCRATCHDIR/contigs_10K.fa \
          --coverage_file $SCRATCHDIR/coverage_table.tsv \
          -b $SCRATCHDIR/concoct_output/

        merge_cutup_clustering.py \
          $SCRATCHDIR/concoct_output/clustering_gt1000.csv \
          > $SCRATCHDIR/concoct_output/clustering_merged.csv

        mkdir -p $SCRATCHDIR/concoct_output/fasta_bins
        extract_fasta_bins.py \
          {input.contigs} \
          $SCRATCHDIR/concoct_output/clustering_merged.csv \
          --output_path $SCRATCHDIR/concoct_output/fasta_bins

        rsync -a $SCRATCHDIR/concoct_output/fasta_bins/ $OUTDIR

        touch {output}
        rm -rf $SCRATCHDIR
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
        "envs/binning_vamb.yaml"
    shell:
        """
        minlen=200000

        ncontigs=$(awk -v minlen="$minlen" '
          /^>/ {{
            if (seqlen >= minlen) n++
            seqlen=0
            next
          }}
          {{ seqlen += length($0) }}
          END {{
            if (seqlen >= minlen) n++
            print n+0
          }}
        ' {input.contigs})
        if [ "$ncontigs" -lt 500 ]; then
            echo "Too few contigs for VAMB, skipping"
            mkdir -p $(dirname {output})
            touch {output}
            exit 0
        fi

        outdir={RESULTS_DIR}/bins/coassembly/vamb
        rm -rf "$outdir"

        vamb \
          --outdir "$outdir" \
          --fasta {input.contigs} \
          --jgi {input.depth} \
          --minfasta {VAMB_MINFASTA} \
          -m {VAMB_MINCONTIGS}

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
        "envs/binning_vamb.yaml"
    shell:
        """
        minlen=200000

        ncontigs=$(awk -v minlen="$minlen" '
          /^>/ {{
            if (seqlen >= minlen) n++
            seqlen=0
            next
          }}
          {{ seqlen += length($0) }}
          END {{
            if (seqlen >= minlen) n++
            print n+0
          }}
        ' {input.contigs})
        if [ "$ncontigs" -lt 500 ]; then
            echo "Too few contigs for VAMB, skipping"
            mkdir -p $(dirname {output})
            touch {output}
            exit 0
        fi

        outdir={RESULTS_DIR}/bins/single/vamb/{wildcards.assembly_type}/{wildcards.sample}
        rm -rf "$outdir"

        vamb \
          --outdir "$outdir" \
          --fasta {input.contigs} \
          --jgi {input.depth} \
          --minfasta {VAMB_MINFASTA} \
          -m {VAMB_MINCONTIGS}

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
        "envs/binning_vamb.yaml"
    shell:
        """
        minlen=200000

        ncontigs=$(awk -v minlen="$minlen" '
          /^>/ {{
            if (seqlen >= minlen) n++
            seqlen=0
            next
          }}
          {{ seqlen += length($0) }}
          END {{
            if (seqlen >= minlen) n++
            print n+0
          }}
        ' {input.contigs})

        if [ "$ncontigs" -lt 500 ]; then
            echo "Too few contigs for VAMB, skipping"
            mkdir -p $(dirname {output})
            touch {output}
            exit 0
        fi


        outdir={RESULTS_DIR}/bins/multi/vamb/{wildcards.assembly_type}/{wildcards.sample}
        rm -rf "$outdir"

        vamb \
          --outdir "$outdir" \
          --fasta {input.contigs} \
          --jgi {input.depth} \
          --minfasta {VAMB_MINFASTA} \
          -m {VAMB_MINCONTIGS}

        touch {output}
        """

############################################
# 5.8. MetaDecoder
############################################

rule metadecoder_coassembly:
    input:
        contigs=f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta",
        bams=expand(
            f"{SCRATCH_MAP}/coassembly/short/{{sample}}.sorted.bam",
            sample=SAMPLES
        )
    params:
        sams=expand(
            f"{SCRATCH_MAP}/coassembly/short/{{sample}}.sorted.sam",
            sample=SAMPLES
        )
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/metadecoder/bins.done")
    threads: config["threads"]
    conda: "envs/binning_metadecoder.yaml"
    run:
        # Generate BAM -> SAM conversion commands in Python
        bam_to_sam_cmds = " && ".join([f"samtools view -h {bam} > {sam}" 
                                       for bam, sam in zip(input.bams, params.sams)])
        
        outdir = f"{RESULTS_DIR}/bins/coassembly/metadecoder"
        
        shell(f"""
        mkdir -p "{outdir}"

        # Convert BAM -> SAM
        {bam_to_sam_cmds}

        # Run MetaDecoder
        metadecoder coverage --threads {threads} -s {params.sams} -o "{outdir}/METADECODER_gsa.COVERAGE"

        metadecoder seed --threads {threads} -f {input.contigs} -o "{outdir}/METADECODER_gsa.SEED"

        metadecoder cluster -f {input.contigs} -c "{outdir}/METADECODER_gsa.COVERAGE" -s "{outdir}/METADECODER_gsa.SEED" -o "{outdir}/METADECODER_coassembly"

        touch {output}
        """)


rule metadecoder_single:
    input:
        contigs=f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bam=f"{SCRATCH_MAP}/single/{{assembly_type}}/{{sample}}.sorted.bam"
    params:
        sam=f"{SCRATCH_MAP}/single/{{assembly_type}}/{{sample}}.sorted.sam"
    output:
        touch(f"{RESULTS_DIR}/bins/single/metadecoder/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda: "envs/binning_metadecoder.yaml"
    run:
        outdir = f"{RESULTS_DIR}/bins/single/metadecoder/{wildcards.assembly_type}/{wildcards.sample}"
        
        shell(f"""
        mkdir -p "{outdir}"

        # Convert BAM -> SAM
        samtools view -h {input.bam} > {params.sam}

        # Run MetaDecoder
        metadecoder coverage --threads {threads} -s {params.sam} -o "{outdir}/METADECODER_gsa.COVERAGE"

        metadecoder seed --threads {threads} -f {input.contigs} -o "{outdir}/METADECODER_gsa.SEED"

        metadecoder cluster -f {input.contigs} -c "{outdir}/METADECODER_gsa.COVERAGE" -s "{outdir}/METADECODER_gsa.SEED" -o "{outdir}/METADECODER_{wildcards.sample}"

        touch {output}
        """)


rule metadecoder_multi:
    input:
        contigs=f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bams=lambda wc: expand(
            f"{SCRATCH_MAP}/multi/{wc.assembly_type}/{wc.sample}/{{rep}}.sorted.bam",
            rep=SAMPLES
        )
    params:
        sams=lambda wc: expand(
            f"{SCRATCH_MAP}/multi/{wc.assembly_type}/{wc.sample}/{{rep}}.sorted.sam",
            rep=SAMPLES
        )
    output:
        touch(f"{RESULTS_DIR}/bins/multi/metadecoder/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda: "envs/binning_metadecoder.yaml"
    run:
        outdir = f"{RESULTS_DIR}/bins/multi/metadecoder/{wildcards.assembly_type}/{wildcards.sample}"
        shell(f"mkdir -p {outdir}")

        # Convert BAM -> SAM for all replicates
        for bam, sam in zip(input.bams, params.sams):
            shell(f"samtools view -h {bam} > {sam}")

        # Now join the SAM paths into a space-separated string for MetaDecoder
        sams_str = " ".join(params.sams)

        # Run MetaDecoder
        shell(f"""
            metadecoder coverage --threads {threads} -s {sams_str} -o {outdir}/METADECODER_gsa.COVERAGE
            metadecoder seed --threads {threads} -f {input.contigs} -o {outdir}/METADECODER_gsa.SEED
            metadecoder cluster -f {input.contigs} -c {outdir}/METADECODER_gsa.COVERAGE -s {outdir}/METADECODER_gsa.SEED -o {outdir}/METADECODER_{wildcards.sample}
            touch {output}
        """)




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
        done = f"{RESULTS_DIR}/bins/coassembly/metabinner/bins.done"
    threads: config["threads"]
    conda:
        "envs/binning_metabinner.yaml"
    shell:
        r"""
        set -euo pipefail

        outdir="{RESULTS_DIR}/bins/coassembly/metabinner"
        mkdir -p "$outdir"

        # Generate k-mer profile
        gen_kmer.py \
          {input.contigs} \
          {config[metabinner][min_contig_len]} \
          {config[metabinner][kmer_size]} \
          "$outdir/kmer.tsv"

        # Run MetaBinner
        run_metabinner.sh \
          -t {threads} \
          -a {input.contigs} \
          -o "$outdir/output" \
          -d {input.coverage} \
          -k "$outdir/kmer.tsv"

        touch {output.done}
        """

rule metabinner_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        coverage = f"{RESULTS_DIR}/depth/single/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        done = f"{RESULTS_DIR}/bins/single/metabinner/{{assembly_type}}/{{sample}}/bins.done"
    threads: config["threads"]
    conda:
        "envs/binning_metabinner.yaml"
    shell:
        r"""
        set -euo pipefail

        outdir="{RESULTS_DIR}/bins/single/metabinner/{wildcards.assembly_type}/{wildcards.sample}"
        mkdir -p "$outdir"

        # Generate k-mer profile
        gen_kmer.py \
          {input.contigs} \
          {config[metabinner][min_contig_len]} \
          {config[metabinner][kmer_size]} \
          "$outdir/kmer.tsv"

        # Run MetaBinner
        run_metabinner.sh \
          -t {threads} \
          -a {input.contigs} \
          -o "$outdir/output" \
          -d {input.coverage} \
          -k "$outdir/kmer.tsv"

        touch {output.done}
        """

rule metabinner_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        coverage = f"{RESULTS_DIR}/depth/multi/{{assembly_type}}/{{sample}}/depth.txt"
    output:
        done = f"{RESULTS_DIR}/bins/multi/metabinner/{{assembly_type}}/{{sample}}/bins.done"
    threads: config["threads"]
    conda:
        "envs/binning_metabinner.yaml"
    shell:
        r"""
        set -euo pipefail

        outdir="{RESULTS_DIR}/bins/multi/metabinner/{wildcards.assembly_type}/{wildcards.sample}"
        mkdir -p "$outdir"

        contig_dir=$(dirname {input.contigs})
        contig_base=$(basename {input.contigs} .fasta)

        # Run gen_kmer in contig directory (MetaBinner hardcodes output path)
        (
          cd "$contig_dir"
          gen_kmer.py \
            "$contig_base.fasta" \
            {config[metabinner][min_contig_len]} \
            {config[metabinner][kmer_size]}
        )

        # Move generated kmer file to outdir
        mv \
          "$contig_dir/${{contig_base}}_kmer_{config[metabinner][kmer_size]}_f{config[metabinner][min_contig_len]}.csv" \
          "$outdir/kmer.tsv"

        # Run MetaBinner
        run_metabinner.sh \
          -t {threads} \
          -a {input.contigs} \
          -o "$outdir/output" \
          -d {input.coverage} \
          -k "$outdir/kmer.tsv"

        touch {output.done}
        """

############################################
# 5.11. SemiBin2
############################################

rule semibin2_coassembly:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta",
        bams = expand(f"{SCRATCH_MAP}/coassembly/short/{{sample}}.sorted.bam", sample=SAMPLES)
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/semibin2/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    params:
        bam_list=lambda wc, input: ' '.join(input.bams)
    shell:
        r"""
        SCRATCH_OUT="{SCRATCH}/semibin2_coassembly_out"
        OUT_FINAL="{RESULTS_DIR}/bins/coassembly/semibin2"

        rm -rf "$SCRATCH_OUT"
        mkdir -p "$SCRATCH_OUT"

        SemiBin2 single_easy_bin \
          -t {threads} \
          -i {input.contigs} \
          -b {params.bam_list} \
          -o "$SCRATCH_OUT" \
          --compression none

        rm -rf "$OUT_FINAL"
        mkdir -p "$OUT_FINAL"
        rsync -a "$SCRATCH_OUT/" "$OUT_FINAL/"

        touch {output}
        """

rule semibin2_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bam = f"{SCRATCH_MAP}/single/{{assembly_type}}/{{sample}}.sorted.bam"
    output:
        touch(f"{RESULTS_DIR}/bins/single/semibin2/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    shell:
        r"""
        SCRATCH_OUT="{SCRATCH}/semibin2_single_{wildcards.assembly_type}_{wildcards.sample}"
        OUT_FINAL="{RESULTS_DIR}/bins/single/semibin2/{wildcards.assembly_type}/{wildcards.sample}"

        rm -rf "$SCRATCH_OUT"
        mkdir -p "$SCRATCH_OUT"

        SemiBin2 single_easy_bin \
          -t {threads} \
          -i {input.contigs} \
          -b {input.bam} \
          -o "$SCRATCH_OUT" \
          --compression none

        rm -rf "$OUT_FINAL"
        mkdir -p "$OUT_FINAL"
        rsync -a "$SCRATCH_OUT/" "$OUT_FINAL/"

        touch {output}
        """

rule semibin2_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bams = lambda wc: expand(
            f"{SCRATCH_MAP}/multi/{wc.assembly_type}/{wc.sample}/{{other}}.sorted.bam",
            other=SAMPLES
        )    
    output:
        touch(f"{RESULTS_DIR}/bins/multi/semibin2/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning.yaml"
    params:
        bam_list=lambda wc, input: ' '.join(input.bams)
    shell:
        r"""
        SCRATCH_OUT="{SCRATCH}/semibin2_multi_{wildcards.assembly_type}_{wildcards.sample}"
        OUT_FINAL="{RESULTS_DIR}/bins/multi/semibin2/{wildcards.assembly_type}/{wildcards.sample}"

        rm -rf "$SCRATCH_OUT"
        mkdir -p "$SCRATCH_OUT"

        SemiBin2 single_easy_bin \
          -t {threads} \
          -i {input.contigs} \
          -b {params.bam_list} \
          -o "$SCRATCH_OUT" \
          --compression none

        rm -rf "$OUT_FINAL"
        mkdir -p "$OUT_FINAL"
        rsync -a "$SCRATCH_OUT/" "$OUT_FINAL/"

        touch {output}
        """

############################################
# 5.12. COMEBin
############################################

rule comebin_coassembly:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/coassembly/short/assembly.fasta",
        bams_dir = directory(f"{SCRATCH_MAP}/coassembly/short")
    output:
        touch(f"{RESULTS_DIR}/bins/coassembly/comebin/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_comebin.yaml"
    shell:
        r"""
        set -euo pipefail

        SCRATCH_OUT="{SCRATCH}/comebin_coassembly_out"
        OUT_FINAL="{RESULTS_DIR}/bins/coassembly/comebin"

        rm -rf "$SCRATCH_OUT"
        mkdir -p "$SCRATCH_OUT"

        run_comebin.sh \
          -t {threads} \
          -a {input.contigs} \
          -o "$SCRATCH_OUT" \
          -p {input.bams_dir}

        rm -rf "$OUT_FINAL"
        mkdir -p "$OUT_FINAL"
        rsync -a "$SCRATCH_OUT/" "$OUT_FINAL/"

        touch {output}
        """

rule comebin_single:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bam = f"{SCRATCH_MAP}/single/{{assembly_type}}/{{sample}}.sorted.bam"
    output:
        touch(f"{RESULTS_DIR}/bins/single/comebin/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_comebin.yaml"
    shell:
        r"""
        set -euo pipefail

        SCRATCH_OUT="{SCRATCH}/comebin_single_{wildcards.assembly_type}_{wildcards.sample}"
        OUT_FINAL="{RESULTS_DIR}/bins/single/comebin/{wildcards.assembly_type}/{wildcards.sample}"
        bamdir={SCRATCH_MAP}/single/{wildcards.assembly_type}/{wildcards.sample}/
        mkdir -p $bamdir
        cp {input.bam} $bamdir

        rm -rf "$SCRATCH_OUT"
        mkdir -p "$SCRATCH_OUT"

        run_comebin.sh \
          -t {threads} \
          -a {input.contigs} \
          -o "$SCRATCH_OUT" \
          -p $bamdir

        rm -rf "$OUT_FINAL"
        mkdir -p "$OUT_FINAL"
        rsync -a "$SCRATCH_OUT/" "$OUT_FINAL/"

        touch {output}
        """

rule comebin_multi:
    input:
        contigs = f"{RESULTS_DIR}/assemblies/single/{{assembly_type}}/{{sample}}/assembly.fasta",
        bams = lambda wc: glob(
            f"{SCRATCH_MAP}/multi/{wc.assembly_type}/{wc.sample}/*.sorted.bam"
        )
    output:
        touch(f"{RESULTS_DIR}/bins/multi/comebin/{{assembly_type}}/{{sample}}/bins.done")
    threads: config["threads"]
    conda:
        "envs/binning_comebin.yaml"
    shell:
        r"""
        set -euo pipefail

        SCRATCH_OUT="{SCRATCH}/comebin_multi_{wildcards.assembly_type}_{wildcards.sample}"
        OUT_FINAL="{RESULTS_DIR}/bins/multi/comebin/{wildcards.assembly_type}/{wildcards.sample}"

        rm -rf "$SCRATCH_OUT"
        mkdir -p "$SCRATCH_OUT"

        run_comebin.sh \
          -t {threads} \
          -a {input.contigs} \
          -o "$SCRATCH_OUT" \
          -p {input.bams}

        rm -rf "$OUT_FINAL"
        mkdir -p "$OUT_FINAL"
        rsync -a "$SCRATCH_OUT/" "$OUT_FINAL/"

        touch {output}
        """

############################################
# 6. Evaluation for rank scoring
############################################

# this BIN_FILES scheme is probably only for Metabat2, other tools have different output. Maybe we have to write a rule for every tool.
import glob
BIN_FILES_SINGLE = lambda wildcards: sorted(
    glob.glob(
        f"{RESULTS_DIR}/bins/single/{wildcards.tool}/{wildcards.assembly_type}/{wildcards.sample}/**/*.fa*",
        recursive=True
    )
)
BIN_FILES_MULTI = lambda wildcards: sorted(glob.glob(f"{RESULTS_DIR}/bins/multi/{wildcards.tool}/{wildcards.assembly_type}/{wildcards.sample}/**/*.fa*", recursive=True))
BIN_FILES_COASSEMBLY = lambda wildcards: sorted(glob.glob(f"{RESULTS_DIR}/bins/coassembly/{wildcards.tool}/**/*.fa*", recursive=True))


############################################
# 6.1 Completeness and Contamination - CheckM 2
############################################

CHECKM2_DB = f"{RESULTS_DIR}/../../database/uniref100.KO.1.dmnd"

rule eval_comp_cont_single:
    input:
        binflag = f"{RESULTS_DIR}/bins/single/{{tool}}/{{assembly_type}}/{{sample}}/bins.done",
        bins = BIN_FILES_SINGLE,
        db = CHECKM2_DB
    output:
        touch(f"{RESULTS_DIR}/eval/single/{{tool}}/{{assembly_type}}/{{sample}}/eval_comp_cont.done")
    log:
        f"logs/{DATASET}/eval/single/{{tool}}/{{assembly_type}}/{{sample}}.comp_cont.log"
    threads: config["threads"]
    conda:
        "envs/evaluation.yaml"
    run:
        # Escape early if input bins are empty
        if not input.bins:
            print(f"[{wildcards.sample}] No bins found for {wildcards.tool} {wildcards.assembly_type}, skipping.")
            shell(f"touch {output[0]}")
        else:
            outdir = f"{RESULTS_DIR}/eval/single/{wildcards.tool}/{wildcards.assembly_type}/{wildcards.sample}_check"
            shell(f"""
                checkm2 predict --threads {threads} --input {' '.join(input.bins)} \
                --output-directory {outdir} --database_path {input.db} --force
                touch {output[0]}
            """)

# Multi-sample Evaluation (short, long, hybrid)
rule eval_comp_cont_multi:
    input:
        binflag = f"{RESULTS_DIR}/bins/multi/{{tool}}/{{assembly_type}}/{{sample}}/bins.done",
        bins = BIN_FILES_MULTI,
        db = CHECKM2_DB
    output:
        touch(f"{RESULTS_DIR}/eval/multi/{{tool}}/{{assembly_type}}/{{sample}}/eval_comp_cont.done")
    log:
        f"logs/{DATASET}/eval/multi/{{tool}}/{{assembly_type}}/{{sample}}.comp_cont.log"
    threads: config["threads"]
    conda:
        "envs/evaluation.yaml"
    run:
        # Escape early if input bins are empty
        if not input.bins:
            print(f"[{wildcards.sample}] No bins found for {wildcards.tool} {wildcards.assembly_type}, skipping.")
            shell(f"touch {output[0]}")
        else:
            outdir = f"{RESULTS_DIR}/eval/multi/{wildcards.tool}/{wildcards.assembly_type}/{wildcards.sample}_check"
            shell(f"""
                checkm2 predict --threads {threads} --input {' '.join(input.bins)} \
                --output-directory {outdir} --database_path {input.db} --force
                touch {output[0]}
            """)

# Co-Assembly Evaluation (short)
rule eval_comp_cont_coassembly:
    input:
        binflag = f"{RESULTS_DIR}/bins/coassembly/{{tool}}/bins.done",
        bins = BIN_FILES_COASSEMBLY,
        db = CHECKM2_DB
    output:
        touch(f"{RESULTS_DIR}/eval/coassembly/{{tool}}/eval_comp_cont.done")
    log:
        f"logs/{DATASET}/eval/coassembly/{{tool}}.comp_cont.log"
    threads: config["threads"]
    conda:
        "envs/evaluation.yaml"
    run:
        # Escape early if input bins are empty
        if not input.bins:
            print(f"No bins found for Coassembly, skipping.")
            shell(f"touch {output[0]}")
        else:
            outdir = f"{RESULTS_DIR}/eval/multi/{wildcards.tool}/{wildcards.assembly_type}/{wildcards.sample}_check"
            shell(f"""
                outdir={RESULTS_DIR}/eval/coassembly/{wildcards.tool}_check
                checkm2 predict --threads {threads} --input {input.bins} --output-directory $outdir --database_path {input.db} --force
                touch {output}
            """)

############################################
# 6.2 tRNA count - Aragorn
############################################

# Single-sample Evaluation (short, long, hybrid)
rule eval_tRNA_single:
    input:
        binflag = f"{RESULTS_DIR}/bins/single/{{tool}}/{{assembly_type}}/{{sample}}/bins.done",
        bins = BIN_FILES_SINGLE
    output:
        final= f"{RESULTS_DIR}/eval/single/{{tool}}/{{assembly_type}}/{{sample}}/eval_tRNA.done"
    log:
        f"logs/{DATASET}/eval/single/{{tool}}/{{assembly_type}}/{{sample}}.tRNA_count.log"
    conda:
        "envs/evaluation.yaml"
    shell:
        """
        eval={RESULTS_DIR}/eval/single/{wildcards.tool}/{wildcards.assembly_type}/{wildcards.sample}
        for filename in {input.bins}; do
                outfile=$eval/tRNA_count.txt
        
                count=$(aragorn -t -w "$filename" | grep -c '^')
                echo -e "$count" >> $outfile
        done
        
        touch {output.final}
        """

# Multi-sample Evaluation (short, long, hybrid)
rule eval_tRNA_multi:
    input:
        binflag = f"{RESULTS_DIR}/bins/multi/{{tool}}/{{assembly_type}}/{{sample}}/bins.done",
        bins = BIN_FILES_MULTI
    output:
        final= f"{RESULTS_DIR}/eval/multi/{{tool}}/{{assembly_type}}/{{sample}}/eval_tRNA.done"
    log:
        f"logs/{DATASET}/eval/multi/{{tool}}/{{assembly_type}}/{{sample}}.tRNA_count.log"
    conda:
        "envs/evaluation.yaml"
    shell:
        """
        eval={RESULTS_DIR}/eval/multi/{wildcards.tool}/{wildcards.assembly_type}/{wildcards.sample}
        for filename in {input.bins}; do
                outfile=$eval/tRNA_count.txt
        
                count=$(aragorn -t -w "$filename" | grep -c '^')
                echo -e "$count" >> $outfile
        done
        
        touch {output.final}
        """

# Co-Assembly Evaluation (short, long, hybrid)
rule eval_tRNA_coassembly:
    input:
        binflag = f"{RESULTS_DIR}/bins/coassembly/{{tool}}/bins.done",
        bins = BIN_FILES_COASSEMBLY
    output:
        final= f"{RESULTS_DIR}/eval/coassembly/{{tool}}/eval_tRNA.done"
    log:
        f"logs/{DATASET}/eval/coassembly/{{tool}}.tRNA_count.log"
    conda:
        "envs/evaluation.yaml"
    shell:
        """
        eval={RESULTS_DIR}/eval/coassembly/{wildcards.tool}
        for filename in {input.bins}; do
                outfile=$eval/tRNA_count.txt
        
                count=$(aragorn -t -w "$filename" | grep -c '^')
                echo -e "$count" >> $outfile
        done
        
        touch {output.final}
        """

############################################
# 6.3 Presence of rRNAs - Barrnap
############################################

# Single-sample Evaluation (short, long, hybrid)
rule eval_rRNA_single:
    input:
        binflag = f"{RESULTS_DIR}/bins/single/{{tool}}/{{assembly_type}}/{{sample}}/bins.done",
        bins = BIN_FILES_SINGLE
    output:
        touch(f"{RESULTS_DIR}/eval/single/{{tool}}/{{assembly_type}}/{{sample}}/eval_rRNA.done")
    log:
        f"logs/{DATASET}/eval/single/{{tool}}/{{assembly_type}}/{{sample}}.rRNA.log"
    threads: config["threads"]
    conda:
        "envs/evaluation.yaml"
    shell:
        """
        for filename in {input.bins}; do
            base=$(basename "$filename" *.fa*)
            outfile="{RESULTS_DIR}/eval/single/{wildcards.tool}/{wildcards.assembly_type}/{wildcards.sample}.${{base}}.barrnap.gff3"   # unique per bin
            barrnap --threads {threads} --quiet "$filename" > "$outfile"

            if grep -E "rRNA.*5S" "$outfile" >/dev/null && \
               grep -E "rRNA.*16S" "$outfile" >/dev/null && \
               grep -E "rRNA.*23S" "$outfile" >/dev/null; then
                all_present=1
            else
                all_present=0
            fi

            echo "$all_present" > "{RESULTS_DIR}/eval/single/{wildcards.tool}/{wildcards.assembly_type}/{wildcards.sample}.${{base}}.all_rRNA_present.txt"
        done
        
        touch {output}
        """

# Multi-sample Evaluation (short, long, hybrid)
rule eval_rRNA_multi:
    input:
        binflag = f"{RESULTS_DIR}/bins/multi/{{tool}}/{{assembly_type}}/{{sample}}/bins.done",
        bins = BIN_FILES_MULTI
    output:
        touch(f"{RESULTS_DIR}/eval/multi/{{tool}}/{{assembly_type}}/{{sample}}/eval_rRNA.done")
    log:
        f"logs/{DATASET}/eval/multi/{{tool}}/{{assembly_type}}/{{sample}}.rRNA.log"
    threads: config["threads"]
    conda:
        "envs/evaluation.yaml"
    shell:
        """
        for filename in {input.bins}; do
            base=$(basename "$filename" *.fa*)
            outfile="{RESULTS_DIR}/eval/multi/{wildcards.tool}/{wildcards.assembly_type}/{wildcards.sample}.${{base}}.barrnap.gff3"   # unique per bin
            barrnap --threads {threads} --quiet "$filename" > "$outfile"

            if grep -E "rRNA.*5S" "$outfile" >/dev/null && \
               grep -E "rRNA.*16S" "$outfile" >/dev/null && \
               grep -E "rRNA.*23S" "$outfile" >/dev/null; then
                all_present=1
            else
                all_present=0
            fi

            echo "$all_present" > "{RESULTS_DIR}/eval/multi/{wildcards.tool}/{wildcards.assembly_type}/{wildcards.sample}.${{base}}.all_rRNA_present.txt"
        done
        
        touch {output}
        """

# Co-Assembly Evaluation (short, long, hybrid)
rule eval_rRNA_coassembly:
    input:
        binflag = f"{RESULTS_DIR}/bins/coassembly/{{tool}}/bins.done",
        bins = BIN_FILES_COASSEMBLY
    output:
        touch(f"{RESULTS_DIR}/eval/coassembly/{{tool}}/eval_rRNA.done")
    log:
        f"logs/{DATASET}/eval/coassembly/{{tool}}.rRNA.log"
    threads: 1
    conda:
        "envs/evaluation.yaml"
    shell:
        """
        for filename in {input.bins}; do
            base=$(basename "$filename" *.fa*)
            outfile="{RESULTS_DIR}/eval/coassembly/{wildcards.tool}.${{base}}.barrnap.gff3"   # unique per bin
            barrnap --threads {threads} --quiet "$filename" > "$outfile"

            if grep -E "rRNA.*5S" "$outfile" >/dev/null && \
               grep -E "rRNA.*16S" "$outfile" >/dev/null && \
               grep -E "rRNA.*23S" "$outfile" >/dev/null; then
                all_present=1
            else
                all_present=0
            fi

            echo "$all_present" > "{RESULTS_DIR}/eval/coassembly/{wildcards.tool}.${{base}}.all_rRNA_present.txt"
        done
        
        touch {output}
        """






