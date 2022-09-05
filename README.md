.###

How to use Snakemake workflow: long-read-chloroplast-assembly
Found at:
https://github.com/aaronphillips7493/long-read-chloroplast-assembly

.###

.###

Use conda:
https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

Make sure you have snakemake and Biopython installed. We use Mamba for increased speed:

a) conda create --name snakemake

b) conda activate snakemake

c) conda install mamba

d) mamba install snakemake biopython

e) snakemake -j 1 --conda-create-envs-only --use-conda (optional - create all conda environments for the pipeline)

.###
Steps
.###

1. Install {The Workflow}:
	git clone https://github.com/aaronphillips7493/long-read-chloroplast-assembly

2. Test {The Workflow}. We provide a test read file containing ONT reads from _Oryza sativa_ ("chloro_assembly/reads/NC_008155.1_single.fasta") and the reference _Oryza sativa_ chloroplast genome ("chloro_assembly/reference/NC_008155.1_single.fasta").
	a) cd long-read-chloroplast-assembly
	b) snakemake --profile profiles/slurm --use-conda --keep-going
	#note: profiles/slurm may not be appropriate for your PC
This test should complete with no errors, and should generate a rotated choloroplast fasta file ("chloro_assembly/{sample}~{assembler}\_chloroplast.fasta") derived from Flye and/or Unicycler. If outside network access is problematic, run downloadReference.sh to download the reference genome. Alternatively, download your reference genome of interest manually and save it into the directory "chloro_assembly/reference/{file_name}\_single.fasta"

3. Run your samples through {The Workflow}.
	a) modify "config.yml":
		i) NCBI_reference_accession = change this to the NCBI accession number for the reference chloroplast genome of interest. Default = NC_008155.1 (O. sativa). If outside network access is problematic, run downloadReference.sh to download the reference genome. Alternatively, download your reference genome of interest manually and save it into the directory "chloro_assembly/reference/{file_name}\_single.fasta"
		ii) my_Email = provide an email address. Neccessary for the automated download of genomes from NCBI - prevents system abuse.
		iii) fast_file = declare whether your read file is in FASTA ("fasta"/"fasta.gz") or FASTQ ("fastq"/"fastq.gz") format.
		iv) randseed = you can supply a seed (integer) here if you wish to re-run read extraction the same way with each redo, or leave this blank to let {The Workflow} randomly generate a seed each time it is run. Change this parameter if assembly fails, and re-run random read subsampling.
		v) numberReads = declare the max. number of reads {The Workflow} will subsample. Change value if assembly fails or to increase or decrease genome coverage. Default = 3000.
		vi) readMinLength = declare the smallest read to be used in the assembly. Default = 5000 bp.
		vii) flyeParameter = tell Flye what kind of reads you are using as input (raw, corrected or HQ ONT long reads, or raw, corrected, or HIFI PacBio long reads). Default = --nano-raw.
		viii) minimap2Parameter = tell minimap2 what kinds of reads you are using (ONT, PacBio, or HIFI). Default = map-ont.
		iX) chloroplastSize = the expected size of the chloroplast genome to be assembled. This is usually set as the size of the reference chloroplast genome. Default = 135000.
		X) cpus = declare the number of CPUs to use. Default = 6.
	b) make sure you have saved your reads in the directory:
		chloro_assembly/reads
	c) make sure you have saved your reference genome in the directory:
		chloro_assembly/reference
	d) run {The Workflow}:
		snakemake --profile profiles/slurm --use-conda --keep-going
