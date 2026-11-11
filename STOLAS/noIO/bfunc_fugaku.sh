#!/bin/bash

# Recompute compaction() from the Nmap files made by bstolas_fugaku.sh,
# recentering on the mu2 peak (functions.cpp).
#
#   Build:  make -f Makefile.fugaku.func
#   Submit: pjsub bfunc_fugaku.sh
#
# Reads : $data/$model/Nfile/Nmap_<NL>_<NF>_<seed>_0_<SIGMA>.bin, $data/$model/logwfile/logw_<...>.dat
# Writes: $data/$model/cmpfile_recenter/compaction_<...>_recenter.dat
#         $data/$model/probabilities_recenter.dat (seed-sorted)
#
# One node, `proc` independent processes (one per CMG), each walking through
# every proc-th seed. probabilities_recenter.dat is sorted by seed at the end.
# parameters.hpp / model.hpp (NLnoise, MODEL, nbias, ...) must match the STOLAS run.

#PJM -L "node=1"
#PJM -L "rscgrp=small"
#PJM -L "elapse=01:00:00"
#PJM -g hp260514
#PJM -x PJM_LLIO_GFSCACHE=/vol0005
#PJM --mpi "proc=4"
#PJM -j
#PJM -S
#PJM --spath "stats/%n.%j.stats"

module load lang/tcsds-1.2.43

# Seeds: ll, ll+1, ..., ll+numdata-1
ll=0
numdata=100
# calPzeta0 (SIGMA) values, same as in bstolas_fugaku.sh
SIGMAS="500."

PROCS=${PJM_MPI_PROC:-1}
THREADS=$((48 / PROCS))
export OMP_NUM_THREADS=$THREADS
export XOS_MMM_L_PAGING_POLICY="demand:demand:demand"

MODEL=$(grep '^#define MODEL' model.hpp | awk '{print $3}')

data=$(grep 'const std::string sdatadir' parameters.hpp \
        | sed 's|//.*||' \
        | sed 's/.*= *"\(.*\)".*/\1/')

if [ "$MODEL" -eq 0 ]; then
    model="chaotic"
elif [ "$MODEL" -eq 1 ]; then
    model="Starobinsky"
elif [ "$MODEL" -eq 2 ]; then
    model="USR"
elif [ "$MODEL" -eq 3 ]; then
    model="hybrid"
else
    echo "Unknown MODEL=$MODEL"
    exit 1
fi

outdir="$data/$model"
logdir="$outdir/log_func"
mkdir -p "$logdir"
mkdir -p "$outdir/cmpfile_recenter"
rm -f "$outdir/probabilities_recenter.dat"

run_worker() {
    local rank=${OMPI_COMM_WORLD_RANK:-${PMIX_RANK:-0}}
    local i
    for ((i = ll + rank; i < ll + numdata; i += PROCS))
    do
    for m in $SIGMAS
    do
        local log="$logdir/${i}_${m%.}.log"
        echo "seed=$i calPzeta=$m rank=$rank threads=$OMP_NUM_THREADS" > "$log"
        SIGMA=$m ./functions $i >> "$log" 2>&1
        local rc=$?
        if [ "$rc" -ne 0 ]; then
            echo "seed=$i calPzeta=$m exit=$rc" >> "$logdir/failed.txt"
        fi
    done
    done
    return 0 # never abort the other ranks on this node
}

export -f run_worker
export ll numdata SIGMAS PROCS logdir

mpiexec -n "$PROCS" bash -c run_worker

# sort the results (columns: seed, mu2hat, k3hat, k_* r_m, zeta_m, Cbar_m, ln W)
sort -n -k1,1 -o "$outdir/probabilities_recenter.dat" "$outdir/probabilities_recenter.dat"
[ -f "$logdir/failed.txt" ] && echo "some seeds failed: see $logdir/failed.txt"
exit 0
