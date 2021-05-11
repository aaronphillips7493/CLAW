###
How to use Snakemake workflow: long-read-chloroplast-assembly
Found at:
https://github.com/aaronphillips7493/long-read-chloroplast-assembly
###

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

3. open script "download_chloro_ref.sh" in a text editor

4. find your chloroplast genome(s) of interest in NCBI and add their NCBI accession numbers to the array in "download_chloro_ref.sh"

e.g. declare -a arr=("NC_041421.1")

If you have multiple target chlorplast genomes (I used 22 rice chloroplast genomes for the rice assembly):

declare -a arr=("KF359913.1"    "KJ830774.1"    "KM881634.1"    "KT992850.1"    "MG383937.1"    "KF359912.1"    "KM881640.1"    "KF359914.1"\
	"KF359915.1"    "KF359918.1"    "KM881641.1"    "JN005831.1"    "KF359921.1"    "KU179220.1"    "KM881643.1"    "KM103375.1"    "KF359911.1"    \
	"KF359919.1"    "JN005832.1"    "KM103369.1"    "AP006728.1"    "AY522329.1")

5. run "download_chloro_ref.sh" using bash in a direcotry that contains chloro_assembly/reference;
this will generate a file containing the chloroplast genome of interest duplicated to mimic the circularisation of a true chloroplast genome, and thus allow long reads to map across the artificial break point introduced into chloroplast assemblies

6. Adjust --genome-size and --threads as necessary under rule raw_flye_chloro_assembly 

7. Adjust the {asm} and {overlap} wild cards under the rule all section of the Snakefile as required; 
if you have three asm values set (e.g. asm=["50","100","300"]) and two overlap values (e.g.overlap=["5000","10000"]) Flye will make 6 (2*3) assemblies.

8. Run snakemake on Snakefile_210510

9. scream into the void as somethign inevitably goes wrong.
