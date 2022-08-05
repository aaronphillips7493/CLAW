.###

How to use Snakemake workflow: long-read-chloroplast-assembly
Found at:
https://github.com/aaronphillips7493/long-read-chloroplast-assembly

.###

.###

Use conda:
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

Make sure you have snakemake and Biopython installed:

a) conda create --name snakemake

b) conda activate snakemake

c) conda install mamba

d) mamba install snakemake biopython

e) snakemake -j 1 --conda-create-envs-only --use-conda (optional - create all conda environments for the pipeline)

.###

#####Still need to fix the yml file situation; try have one yml file for each rule in Snakefile#####

#####Still need to find a better way of selecting PacBio or ONT mapping options for read mapping#####

#####Still need to find a better way to enter sample names as the prefix#####

#####Still need to find a better way to enter genome size#####
  
.###
Steps
.###

1. Add your unique sample prefix(s) to the SAMPLES section of the Snakefile. 

e.g. for me, I had long reads from four populations of O. australiensis, and the fastq read files were named as follows:

Oaustraliensis-300131-flowcell-1-SQK-LSK109_guppy303_all.fastq,

Oaustraliensis-keepriver-flowcell-2-SQK-LSK109_guppy303_all.fastq,

Oaustraliensis_keepriver_1g_SQK-LSK109_guppy303_all.fastq,

etc...

So, to the SAMPLES section, I add:

"Oaustraliensis-300131-flowcell-1-SQK-LSK109",

"Oaustraliensis-keepriver-flowcell-2-SQK-LSK109",

"Oaustraliensis_keepriver_1g_SQK-LSK109",

etc...  

2. The first step of the pipeline converts a guppy base-called fastq file with the suffix *_guppy303_all.fastq to a fastq.gz file. 
Have your ONT long read fastq file saved in a directory called backup/{prefix}/, where prefix is the prefix you saved to the SAMPLES list.  

If you already have your ONT long reads saved as *_guppy303_all.fastq.gz, move them to a directory called working_raw_base_calls 
(need to make pipeline amenable to PacBio reads too; need to make the SAMPLES bit more general [no need to add your own prefix to SAMPLES] and remove the "*_guppy303_all.fastq" thing)

3. find your chloroplast genome(s) of interest in NCBI and add their NCBI accession numbers to NCBI.remote() in the rule "download_chloro_genome"

e.g. NCBI.remote("NC_041421.1", db="nuccore")

###I dont't know if this will work for multiple reference genomes yet###

4. Adjust:
	in rule minimap2_plastid: select ONT or PacBio reads (map-ont for ONT or map-pb for PacBio) (((Need to automate this process))) 
	in rule raw_flye_chloro_assembly: select PacBio as input; set --genome-size and --threads as necessary 

5. Adjust the {asm} and {overlap} wild cards under the rule all section of the Snakefile as required; 
if you have three asm values set (e.g. asm=["50","100","300"]) and two overlap values (e.g.overlap=["5000","10000"]) Flye will make 6 (2*3) assemblies.

6. Run snakemake on Snakefile_210510
e.g. snakemake --use-conda --profile profiles/slurm -s Snakefile_210510

7. scream into the void as somethign inevitably goes wrong.
