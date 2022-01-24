TARGET=${TARGET:-all}

QSUB="qsub -q {cluster.queue} -l ncpus={threads} -l jobfs={cluster.hdd}"
QSUB="$QSUB -l walltime={cluster.time} -l mem={cluster.mem} -N {rule} -l storage={cluster.storage}"
QSUB="$QSUB -l wd -j oe -P {cluster.project}"

snakemake                                                          \
    -j 1000                                                        \
    --max-jobs-per-second 2                                        \
    --cluster-config cluster-configs/pbs.yaml                      \
    --local-cores ${PBS_NCPUS:-1}                                  \
    --js profiles/pbs/jobscript.sh                                 \
    --nolock                                                       \
    --rerun-incomplete                                             \
    --keep-going                                                   \
    --cluster "$QSUB"                                              \
    "$TARGET"
