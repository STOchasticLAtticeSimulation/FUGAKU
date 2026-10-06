#!/bin/bash

# Run STOLAS for many noise seeds ($i) as a Fugaku bulk job.
#
#   Submit (numdata=100, 4 procs/node -> 25 sub-jobs):
#     pjsub --bulk --sparam "0-24" bstolas_fugaku.sh
#   General form: --sparam "0-$(( (numdata + procs - 1) / procs - 1 ))"
#
# Layout: each sub-job takes 1 node and runs `proc` independent STOLAS processes
# side by side (one seed each, no communication), each with 48/proc OpenMP
# threads. proc=4 means one process per CMG (12 cores + 8GB HBM each). Runs
# are independent, so more processes x fewer threads is more efficient than one
# 48-thread process, as long as one run fits in 32GB/proc of memory.
# To change it, edit `--mpi "proc=N"` below (N = 1, 2, 4, 12, 48 ...).
# Sub-jobs are scheduled independently: a failed/timed-out seed doesn't affect the others.
#
# Without --bulk (plain `pjsub bstolas_fugaku.sh`) this single node runs ALL
# seeds, `proc` at a time, one after another -- fine for tests, but for
# numdata=100 at 256^3 it will not fit in `elapse`.

#PJM -L "node=1"
#PJM -L "rscgrp=small"
#PJM -L "elapse=05:00:00"
#PJM -g hp260514
#PJM -x PJM_LLIO_GFSCACHE=/vol0005
#PJM --mpi "proc=4"
#PJM -j
#PJM -S
#PJM --spath "stats/%n.%j.%b.stats"

module load lang/tcsds-1.2.43

# Seeds: ll, ll+1, ..., ll+numdata-1
ll=0
numdata=100
# calPzeta0 (SIGMA) values, run one after another for each seed
SIGMAS="500."

PROCS=${PJM_MPI_PROC:-1}
THREADS=$((48 / PROCS))
if [ -n "$PJM_BULKNUM" ]; then
    BULK=$PJM_BULKNUM
    STRIDE=$numdata # one seed per rank
else
    BULK=0
    STRIDE=$PROCS   # plain job: each rank walks through every PROCS-th seed
    echo "PJM_BULKNUM is not set: running all $numdata seeds on this node (submit with --bulk --sparam to spread them)"
fi

export OMP_NUM_THREADS=$THREADS
# Allocate memory on first touch instead of pre-paging everything at startup
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

logdir="$data/$model/log"
mkdir -p "$data"
mkdir -p "$data/$model"
mkdir -p "$data/$model/animation"
mkdir -p "$logdir"

# Each mpiexec rank is one STOLAS process. mpiexec is only used as a launcher
# that places the ranks on separate CMGs; STOLAS itself doesn't use MPI.
run_worker() {
    local rank=${OMPI_COMM_WORLD_RANK:-${PMIX_RANK:-0}}
    local i
    for ((i = ll + BULK * PROCS + rank; i < ll + numdata; i += STRIDE))
    do
    for m in $SIGMAS
    do
        local log="$logdir/${i}_${m%.}.log"
        echo "seed=$i calPzeta=$m bulk=$BULK rank=$rank threads=$OMP_NUM_THREADS" > "$log"
        taskset -cp $$ >> "$log" 2>&1 # CPU binding, to check the placement
        SIGMA=$m ./STOLAS $i >> "$log" 2>&1
        local rc=$?
        if [ "$rc" -ne 0 ]; then
            echo "seed=$i calPzeta=$m exit=$rc" >> "$logdir/failed.txt"
        fi
    done
    done
    return 0 # never abort the other ranks on this node
}

export -f run_worker
export ll numdata SIGMAS PROCS BULK STRIDE logdir

mpiexec -n "$PROCS" bash -c run_worker
