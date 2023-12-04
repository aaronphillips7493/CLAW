![CLAW](https://github.com/aaronphillips7493/CLAW/blob/main/profiles/CLAW_icon.png?raw=true)
---------------------------------------------------------

About CLAW
---------------------------------------------------------

CLAW (Chloroplast Long-read Assembly Workflow) is an mostly-automated Snakemake-based workflow for the assembly of chloroplast genomes. CLAW uses chloroplast long-reads, which are baited out of larger read libraries (e.g., an Oxford Nanopore Technologies MinION read library derived from photosynthetic tissue), for assembly with Flye and/or Unicycler. CLAW was designed with the novice bioinformatician in mind - it is easy to install and easy to use, requiring only minimal user input.

---------------------------------------------------------

Download and Install
---------------------------------------------------------

1. Clone Git repository:
	
	git clone https://github.com/aaronphillips7493/CLAW

2. Move into the directory:
	
		cd CLAW

3. Install conda (if not already present):

https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html

4. Make sure you have snakemake and Biopython installed. We use Mamba for increased speed:

	a) Create a conda environment:

		conda create --name snakemake

	b) Activate the conda environment:

		conda activate snakemake

	c) Install mamba in the environment:

		conda install mamba

	d) Install snakemake and biopython in the environment:

		mamba install snakemake biopython

	e) OPTIONAL: create all environments needed for CLAW:

		snakemake -j 1 --conda-create-envs-only --use-conda

---------------------------------------------------------

Steps
---------------------------------------------------------

1. Test CLAW. We provide a test read file containing ONT reads from _Oryza sativa_ ("chloro_assembly/reads/DRR196880_subset.fastq") and the reference _Oryza sativa_ chloroplast genome ("chloro_assembly/reference/NC_008155.1_single.fasta").

		snakemake --profile profiles/slurm --use-conda --keep-going
	
	#note: profiles/slurm may not be appropriate for your PC. If this is the case, please run CLAW after specifying your profile (slurm, pbs, or local). 

This test should complete with no errors, and should generate a rotated choloroplast fasta file ("chloro_assembly/{sample}~{assembler}\_chloroplast.fasta") derived from Flye and/or Unicycler. If outside network access is problematic, run downloadReference.sh to download the reference genome. Alternatively, download your reference genome of interest manually and save it into the directory 

	"chloro_assembly/reference/{NCBI_accession_number}_single.fasta"

2. Run your samples through CLAW.

	a) modify "config.yml":
	
		i) ncbi_reference_accession = change this to the NCBI accession number for the reference chloroplast genome of interest. Default = NC_008155.1 (O. sativa). If outside network access is problematic, run downloadReference.sh to download the reference genome. Alternatively, download your reference genome of interest manually and save it into the directory "chloro_assembly/reference/{file_name}\_single.fasta". Make sure to change this parameter even if you do not use CLAW to download the reference genome.
		
		ii) my_Email = provide an email address. Neccessary for the automated download of genomes from NCBI - prevents system abuse.
		
		iii) fast_file = declare whether your read file is in FASTA ("fasta"/"fasta.gz") or FASTQ ("fastq"/"fastq.gz") format.
		
		iv) randseed = you can supply a seed (integer) here if you wish to re-run read extraction the same way with each redo, or leave this blank to let CLAW randomly generate a seed each time it is run. Change this parameter if assembly fails, and re-run random read subsampling.
		
		v) numberReads = declare the max. number of reads CLAW will subsample. Change value if assembly fails or to increase or decrease genome coverage. Default = 3000.
		
		vi) readMinLength = declare the smallest read to be used in the assembly. Default = 5000 bp.
		
		vii) flyeParameter = tell Flye what kind of reads you are using as input (raw, corrected or HQ ONT long reads, or raw, corrected, or HIFI PacBio long reads). Default = --nano-raw.
		
		viii) minimap2Parameter = tell minimap2 what kinds of reads you are using (ONT, PacBio, or HIFI). Default = map-ont.
		
		iX) chloroplastSize = the expected size of the chloroplast genome to be assembled. This is usually set as the size of the reference chloroplast genome. Default = 135000.
		
		X) cpus = declare the number of CPUs to use. Default = 4.
		
	b) make sure you have saved your reads in the directory:
		
		chloro_assembly/reads
		
	c) make sure you have saved your reference genome in the directory:
		
		chloro_assembly/reference
		
	d) run CLAW (remember to specify the correct profile based on our system!):
		
		snakemake --profile profiles/slurm --use-conda --keep-going
		
3. If {The Workflow} fails, try modifying "randSeed" and/or "numberReads" in config.yml. You will need to delete the files in the following directories to re-run {The Workflow}:

		chloro_assembly/assemblies
	
		chloro_assembly/subReads
	
		chloro_assembly/alignments
	
		chloro_assembly/dotPlots

You may also need to delete the file:

	chloro_assembly/{sample}~{assembler}_chloroplast.fasta

---------------------------------------------------------

Notes
---------------------------------------------------------

User can check read mapping and/or coverage of the genomes (whether that be to the single or 'circular' reference genome, or to the CLAW-assembled genome) by loading the genome of interest and its associated BAM/BigWig files into any genome browser (e.g., IGV). Read depth may help you diagnose issues in your assembly, if you have any.

While the read library used for chloroplast genome assembly is enriched for chloroplast reads, there is a possibility that the library will also contain mitochondrial reads because of high sequence similarity between mitochondrial and chloroplast genomes. Thus, it is possible for CLAW to assemble mitochondrial (likely fragmented) contigs too. CLAW makes no attempt to resolve this, so users will have to manually check the assembled contigs for similarity to chloroplast and/or mitochondrial sequences.

While we report on the successful runs of CLAW here, there were several instances in which CLAW failed to complete or when chloroplast genome assembly failed. Typically, CLAW failed to complete due to insufficient memory allocation for certain steps of CLAW. These steps were commonly the read mapping and genome assembly steps. In this case, we were able to overcome this issue by modifying the memory allocation in the “cluster-configs/default.yaml” file. Future releases should aim to integrate dynamic memory allocations to further minimise user input for the successful completion of CLAW. When genome assembly failed, we found that re-running CLAW with a different seed (randomly generated as part of the CLAW, unless specified by the user in “config.yaml”) resulted in successful completion of genome assembly. We also found it necessary to, on occasion, adjust the number of chloroplast reads used for genome assembly as too many reads (and therefore too high genome coverage) could cause the assemblers to fail. Thus, by modifying the random seed and the number of reads used for assembly we were able to generate assemblies for all 19 test datasets. 
