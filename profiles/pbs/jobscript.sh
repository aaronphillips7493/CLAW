#!/bin/bash -l
# properties = {properties}

export TMPDIR=$PBS_JOBFS

eval "$(/g/data/xe2/gadi/conda/bin/conda shell.zsh hook)"
conda activate sf-snake
module load minimap2 samtools
#PATH=$PATH:/g/data/xe2/scott/gadi_modules/MUMmer3.23

set -ueo pipefail
set -x

{exec_job}
