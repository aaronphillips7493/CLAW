import json
configfile: "config.yml"

shell.executable("/bin/bash")
shell.prefix("set -euo pipefail; ")

##########
# local rules
# these rules will be performed within the main job when using a cluster
localrules: all,
    rotate_chloroplast,
    dot_plot,
    sub_sample,
    extract_aligned_reads,
    index_reference,
    double_chloro_genome,
    rotate_chloroplast

# need internet access
#   download_chloro_genome

# big jobs, when using HPC run in queue.
#   assemble
#   align

##########
# Generate random seed for subsampling reads
import random
RAND_SEED= config["rand_seed"]
try: RAND_SEED
except NameError: RAND_SEED = random.seed()

##########
# avoid recursion loop
wildcard_constraints:
    sample="[^/]+",

##########
# Find all fasta/fastq/fasta.gz/fastq.gz files
# continaing reads that we want to build our genome from
# and attempt assembly of chloroplast genome
#
samples = glob_wildcards("chloro_assembly/reads/{file}." + config["fast_file"]).file
#print(samples)

rule all:
    input:
        "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_index.fasta",
        expand("chloro_assembly/dotPlots/{sample}_{assembler}.png", sample=samples, assembler=config["assembler"],),
        expand("chloro_assembly/alignments/{sample}_subset_vs_{assembler}_assembly_sorted.bw", sample=samples, assembler=config["assembler"],),
        expand("chloro_assembly/alignments/{sample}_subset_vs_ref_sorted.bw", sample=samples,),
        #expand("chloro_assembly/{sample}~chloroplast.fasta", sample=samples)

##########
# We have a chloroplast genome.
# Now we want to "rotate" our assemblies such that the linear sequence
# of our genome starts at the same position as the used reference genome
# And as our assemblies may have multiple contigs we need to rotate all.
#
rule rotate_chloroplast:
    input:
        index = "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_index.fasta",
        genome = "chloro_assembly/assemblies/{sample}_{assembler}/assembly.fasta"
    output:
        "chloro_assembly/{sample}~{assembler}_chloroplast.fasta",
    params:
        dir = "chloro_assembly/{sample}~{assembler}_tmp",
        name = "{sample}~{assembler}_chloroplast"
    conda:
        "envs/rotate_chloroplast.yml"
    benchmark:
        "chloro_assembly/benchmark/rotate_chloroplast/{sample}_{assembler}_benchmark.txt"
    shell:
        """
        if [ -s {input.genome} ]
        then
            mkdir -p {params.dir}
            for contig in `grep '^>' {input.genome} | sed -e 's/>//g'`
            do
                echo $contig > {params.dir}/tmp
                seqtk subseq {input.genome} {params.dir}/tmp > {params.dir}/$contig.fasta

                nucmer --maxmatch {input.index} {params.dir}/$contig.fasta -p {params.dir}/out
                show-coords -THrd {params.dir}/out.delta > {params.dir}/out.coords
                start=`sort -k6,6hr {params.dir}/out.coords | head -n 1| cut -f3`
                echo ">$contig" >> {output}
                echo "$start XXX"
                if [ $start == 1 ]
                then
                    grep -v '^>' {params.dir}/$contig.fasta | tr -d '\n' >> {output}
                    echo "" >> {output} 
                elif [ ! -z $start ]
                then
                    grep -v '^>' {params.dir}/$contig.fasta | tr -d '\n' > {params.dir}/temp.fasta
                    cut -c ${{start}}- {params.dir}/temp.fasta > {params.dir}/start.fasta
                    cut -c -$[start-1] {params.dir}/temp.fasta > {params.dir}/end.fasta
                    cat {params.dir}/start.fasta {params.dir}/end.fasta | tr -d '\n' >> {output}
                    echo "" >> {output}
                else
                    grep -v '^>' {params.dir}/$contig.fasta | tr -d '\n' >> {output}
                    echo "" >> {output}
                fi
                rm -rf {params.dir}/*
            done
            rm -rf {params.dir}
        else
            touch {output}
        fi
        """


##########
# After completino of assembly we want to align the new chloroplast to
# the refernec egnome and generate a dot plot.
# The user then needs to view and decide if the assembly contains the chloroplast.
# In the case of a assembly having >1 contig this will allow the user to select the choloropast genome
# i.e. the one that best alignes and is the correct size.

# add to conda config - conda install -c bioconda gnuplot
rule dot_plot:
    input:
        reference = "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_single.fasta",
        query = "chloro_assembly/{sample}~{assembler}_chloroplast.fasta"
    output:
        "chloro_assembly/dotPlots/{sample}_{assembler}.png"
    params:
        "chloro_assembly/dotPlots/{sample}_{assembler}"
    conda:
        "envs/dot_plot.yml"
    benchmark:
        "chloro_assembly/benchmark/dot_plot/{sample}_{assembler}_benchmark.txt"
    shell:
        """
        nucmer {input.reference} {input.query} -p {params}
        # delta-filter -1 -i 50 {params}.delta > {params}.1delta
        mummerplot --fat --png --large {params}.delta -p {params}
        rm {params}.filter {params}.fplot {params}.gp {params}.rplot
        """

##########
# To see which parts of the assembled genomes are covered by the read subset
# View in a genome browser (e.g. IGV; load assembled genome + .bw)

rule bam_coverage_assembly:
    input:
        "chloro_assembly/alignments/{sample}_subset_vs_{assembler}_assembly_sorted.bam",
    output:
        bigwig = "chloro_assembly/alignments/{sample}_subset_vs_{assembler}_assembly_sorted.bw"
    conda:
        "envs/align.yml"
    benchmark:
        "chloro_assembly/benchmark/bam_coverage_assembly/{sample}_{assembler}_benchmark.txt"
    shell:
        """
        bamCoverage -b {input} -o {output.bigwig}
        """

rule sort_assembly_alignment:
    input:
        "chloro_assembly/alignments/{sample}_subset_vs_{assembler}_assembly.bam",
    output:
        bam = "chloro_assembly/alignments/{sample}_subset_vs_{assembler}_assembly_sorted.bam",
    conda:
        "envs/align.yml"
    benchmark:
        "chloro_assembly/benchmark/sort_assembly_alignment/{sample}_{assembler}_benchmark.txt"
    shell:
        """
        samtools sort -o {output.bam} {input}
        samtools index {output.bam}
        """

rule align_subset_to_assembly:
    input:
        fastFile = "chloro_assembly/subReads/{sample}~assemble.fasta",
        assembly = "chloro_assembly/assemblies/{sample}_{assembler}/assembly.fasta"
    output:
        bam = temp("chloro_assembly/alignments/{sample}_subset_vs_{assembler}_assembly.bam"),
    threads:
        config["cpus"]
    conda:
        "envs/align.yml"
    benchmark:
        "chloro_assembly/benchmark/align_subset_to_assembly/{sample}_{assembler}_benchmark.txt"
    shell:
        """
        minimap2 -ax {config[minimap2_parameter]} -t {threads} {input.assembly} {input.fastFile} \
            | samtools view -b -F 4 -@ {threads} \
            > {output.bam}
        """

##########
# Assemble our subset of chloroplast reads.
# If sucessful this will produce our chloroplast genome
# If not we will have to rerun using a different random seed and/or different number of reads.
#
rule unicycler_assemble:
    input:
        "chloro_assembly/subReads/{sample}~assemble.fasta"
    params:
        assembly = "chloro_assembly/assemblies/{sample}_unicycler",
    output:
        "chloro_assembly/assemblies/{sample}_unicycler/assembly.fasta",
    threads:
        config["cpus"]
    conda:
        "envs/assemble.yml"
    benchmark:
        "chloro_assembly/benchmark/unicycler_assemble/{sample}_unicycler_benchmark.txt"
    shell:
        """
        mkdir -p chloro_assembly/assemblies/
        unicycler -l {input} -o {params.assembly}
        """

rule flye_assemble:
    input:
        "chloro_assembly/subReads/{sample}~assemble.fasta"
    params:
        flye_assembly = "chloro_assembly/assemblies/{sample}_flye",
    output:
        "chloro_assembly/assemblies/{sample}_flye/assembly.fasta",
    threads:
        config["cpus"]
    conda:
        "envs/assemble.yml"
    benchmark:
        "chloro_assembly/benchmark/flye_assemble/{sample}_flye_benchmark.txt"
    shell:
        """
        mkdir -p chloro_assembly/assemblies/
        flye --threads {threads} --genome-size {config[chloroplast_size]} -o {params.flye_assembly} {config[flye_parameter]} {input}
        """

##########
# To see which parts of the ref genome are covered by the read subset
# View in a genome browser (e.g. IGV; load ref genome + .bw)

rule bam_coverage_ref:
    input:
        "chloro_assembly/alignments/{sample}_subset_vs_ref_sorted.bam",
    output:
        bigwig = "chloro_assembly/alignments/{sample}_subset_vs_ref_sorted.bw"
    conda:
        "envs/align.yml"
    benchmark:
        "chloro_assembly/benchmark/bam_coverage_ref/{sample}_benchmark.txt"
    shell:
        """
        bamCoverage -b {input} -o {output.bigwig}
        """

rule sort_ref_alignment:
    input:
        "chloro_assembly/alignments/{sample}_subset_vs_ref.bam",
    output:
        bam = "chloro_assembly/alignments/{sample}_subset_vs_ref_sorted.bam",
    conda:
        "envs/align.yml"
    benchmark:
        "chloro_assembly/benchmark/sort_ref_alignment/{sample}_benchmark.txt"
    shell:
        """
        samtools sort -o {output.bam} {input}
        samtools index {output.bam}
        """

rule align_subset_to_ref:
    input:
        fastFile = "chloro_assembly/subReads/{sample}~assemble.fasta",
        reference = "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_circular.fasta"
    output:
        bam = "chloro_assembly/alignments/{sample}_subset_vs_ref.bam",
    threads:
        config["cpus"]
    conda:
        "envs/align.yml"
    benchmark:
        "chloro_assembly/benchmark/align_subset_to_ref/{sample}_benchmark.txt"
    shell:
        """
        minimap2 -ax {config[minimap2_parameter]} -t {threads} {input.reference} {input.fastFile} \
            | samtools view -b -F 4 -@ {threads} \
            > {output.bam}
        """

##########
# From our chloroplast read set randomly sample X reads.
# These will be used to assemble our chloroplast genome.
# We sample to reduce coverage, speeding up assembly.
#
rule sub_sample:
    input:
        "chloro_assembly/subReads/{sample}~all.fasta"
    output:
        fastFile = "chloro_assembly/subReads/{sample}~assemble.fasta"
    params:
        random_seed = RAND_SEED,
    conda:
        "envs/sub_sample.yml"
    benchmark:
        "chloro_assembly/benchmark/sub_sample/{sample}_benchmark.txt"
    shell:
        """
        echo {params.random_seed}
        seqtk sample -s {params.random_seed} {input} {config[number_reads]} > {output}
        """

##########
# After alignement we want to ID and extract chloroplast reads,
# ie. reads that aligned to the reference chloroplast.
#
rule extract_aligned_reads:
    input:
        bam = "chloro_assembly/alignments/{sample}.bam",
        fastFile = "chloro_assembly/reads/{sample}." + config["fast_file"]
    output:
        list = "chloro_assembly/alignments/{sample}~all.lst",
        fastFile = "chloro_assembly/subReads/{sample}~all.fasta"
    conda:
        "envs/extract_aligned_reads.yml"
    benchmark:
        "chloro_assembly/benchmark/extract_aligned_reads/{sample}_benchmark.txt"
    shell:
        """
        samtools view {input.bam} | cut -f1 | sort | uniq > {output.list}
        seqtk subseq {input.fastFile} {output.list} \
            | bioawk -c fastx \
                'length($seq) > {config[read_min_length]} && length($seq) < {config[chloroplast_size]} \
                {{print \">\"$name\"\\n\"$seq}}' > {output.fastFile}
        """

##########
# Align reads to reference genome.
# In order to find chloroplast reads.
#
rule align:
    input:
        fastFile = "chloro_assembly/reads/{sample}." + config["fast_file"],
        reference = "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_circular.fasta"
    output:
        bam = "chloro_assembly/alignments/{sample}.bam",
    threads:
        config["cpus"]
    conda:
        "envs/align.yml"
    benchmark:
        "chloro_assembly/benchmark/align/{sample}_benchmark.txt"
    shell:
        """
        minimap2 -ax {config[minimap2_parameter]} -t {threads} {input.reference} {input.fastFile} \
            | samtools view -b -F 4 -@ {threads} \
            | samtools sort -o {output.bam}
        """

##########
# We want to have our genome start at the same point.
# Get the first 10 Kbp of the reference genome, use this
# to roate our new chloroplast after assembly.
#
rule index_reference:
    input:
        "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_circular.fasta"
    output:
        "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_index.fasta"
    conda:
        "envs/index_reference.yml"
    benchmark:
        "chloro_assembly/benchmark/index_reference/"+config["NCBI_reference_accession"]+"_benchmark.txt"
    shell:
        """
        awk 'NR == 1 {{print substr($1,2,length($1)), \"0\", \"10000\"}}' {input} > chloro_assembly/reference/index.bed
        seqtk subseq {input} chloro_assembly/reference/index.bed > {output}
        rm chloro_assembly/reference/index.bed
        """

##########
# "circularise" our reference genome.
#  Genomes are stored as linear fasta files. Chloroplasts are circular.
#  To better allow reads to align to the break in our circular genome we double up the genome
#
rule double_chloro_genome:
    input:
        "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_single.fasta"
    output:
        "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_circular.fasta"
    benchmark:
        "chloro_assembly/benchmark/double_chloro_genome/"+config["NCBI_reference_accession"]+"_benchmark.txt"
    shell:
        """
        head -n 1 {input} > {output}
        tail -n +2 {input} | tr -d '\n' >> {output}
        tail -n +2 {input} | tr -d '\n' >> {output}
        """


##########
# Download reference chloroplast genome
#  This genome is used to find chloroplast genomes within your read set
#  Need to install BioPython before this rule runs --> installed as part of Snakemake env.
#  If the computer or remote HPC you are working from has restricted access to external
#  data sources, you will need to comment out each line of this rule and manually download
#  each reference genome of interest and save it so it is identical to 'output' below

from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
NCBI = NCBIRemoteProvider(email=config["my_email"]) # email required by NCBI to prevent abuse
rule download_chloro_genome:
    input:
        NCBI.remote(config["NCBI_reference_accession"] +".fasta", db="nuccore")
    output:
        "chloro_assembly/reference/" + config["NCBI_reference_accession"] + "_single.fasta"
    benchmark:
        "chloro_assembly/benchmark/download_chloro_genome/" + config["NCBI_reference_accession"] + "_benchmark.txt"
    shell:
        """
        mv {input} {output}
        """




##########
##########
# Old rules - working but not in use
##########
##########

###rule rotate_chloroplast:
###    input:
###        index = "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_index.fasta",
###        genome = "chloro_assembly/validated/{sample}~chloroplast.fasta"
###    output:
###        "chloro_assembly/{sample}~chloroplast.fasta"
###    params:
###        dir = "chloro_assembly/{sample}~tmp",
###        name = "{sample}~chloroplast"
###    conda:
###        "envs/rotate_chloroplast.yml"
###    benchmark:
###        "chloro_assembly/benchmark/rotate_chloroplast/{sample}_benchmark.txt"
###    shell:
###        """
###        if [ -s {input.genome} ]
###        then
###            mkdir {params.dir}
###            nucmer --maxmatch {input.index} {input.genome} -p {params.dir}/out
###            show-coords -THrd {params.dir}/out.delta > {params.dir}/out.coords
###            start=`sort -k6,6hr {params.dir}/out.coords | head -n 1| cut -f3`
###            echo ">{params.name}" > {output}
###            grep -v '^>' {input.genome} | tr -d '\n' > {params.dir}/temp.fasta
###            cut -c ${{start}}- {params.dir}/temp.fasta > {params.dir}/start.fasta
###            cut -c -$[start-1] {params.dir}/temp.fasta > {params.dir}/end.fasta
###            cat {params.dir}/start.fasta {params.dir}/end.fasta | tr -d '\n' >> {output}
###            rm -rf {params.dir}
###        else
###            touch {output}
###        fi
###        """

##########
# After the contig size has been found to be correct, we want to
# align our potential new chloroplast genome to the reference
# and check it's percent identity. If the identity is within the given range
# we have our new chloroplast
#
###rule check_identity:
###    input:
###        query = "chloro_assembly/checkSize/{sample}.fasta",
###        reference = "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_single.fasta"
###    params:
###        check = "chloro_assembly/checkIdentity/{sample}"
###    output:
###        "chloro_assembly/validated/{sample}~chloroplast.fasta"
###    conda:
###        "envs/check_identity.yml"
###    benchmark:
###        "chloro_assembly/benchmark/check_identity/{sample}_benchmark.txt"
###    shell:
###        """
###        echo {input.query}
###        mkdir -p chloro_assembly/checkIdentity
###        if [ -s {input.query} ]
###        then
###            dnadiff {input.reference} {input.query} -p {params.check}
###            rm {params.check}*.1coords {params.check}*.1delta {params.check}*.delta {params.check}*.mcoords {params.check}*.mdelta {params.check}*.qdiff {params.check}*.rdiff {params.check}*.snps
###            grep 'AvgIdentity' {params.check}.report | awk 'NR == 1{{if(($2/100) >= config[chloroplastIdentityPC]){{print $2}}}}' > {params.check}.good
###            [ -s {params.check}.good ] && cp {input.query} {output}
###        else
###            touch {output}
###        fi
###        """

##########
# After assembly we want to check all assebmled contigs to see
# if their size is approximetely correct.
# If we dont get a contig of the correct size range the prodcued fasta will be empty.
#
###rule check_size:
###    input:
###        reference = "chloro_assembly/reference/"+config["NCBI_reference_accession"]+"_single.fasta",
###        query = "chloro_assembly/assemblies/{sample}/assembly.fasta"
###    output:
###        "chloro_assembly/checkSize/{sample}.fasta"
###    conda:
###        "envs/check_size.yml"
###    benchmark:
###        "chloro_assembly/benchmark/check_size/{sample}_benchmark.txt"
###    shell:
###        """
###        bioawk -c fastx 'BEGIN{{lower = {config[chloroplastSize]} * {config[chloroplastSizePC]}; upper = {config[chloroplastSize]} * (2-{config[chloroplastSizePC]})}}{{if(length($seq) > lower && length($seq) < upper){{print \">\"$name\"\\n\"$seq}}}}' {input.query} > {output}
###        """
