#!/bin/bash

# ---------------------------------------------------------------------------
# pjsub batch script for Fugaku (Fujitsu A64FX, 1 node / 48 cores, OpenMP only
# -- no MPI, so a single node is all this job can use).
#
# Submit with:
#   pjsub bstolas_fugaku.sh
# ---------------------------------------------------------------------------

#PJM -L "node=1"
#PJM -L "rscgrp=int"
#PJM -L "elapse=05:00:00"
#PJM -g hp260514
#PJM -x PJM_LLIO_GFSCACHE=/vol0005
#PJM -j
#PJM -S
#PJM -o job_out

# Match the module used to build STOLAS with Makefile.fugaku.
module load lang/tcsds-1.2.43

# A64FX has 48 usable compute cores per node (assistant cores are reserved by
# the OS and excluded automatically).
export OMP_NUM_THREADS=48

# For the FFTW thread-count sweep diagnostic, use bstolas_fugaku_sweep.sh
# instead (a separate script, not a flag/env-var on this one): pjsub
# --interact doesn't reliably forward the submitting shell's environment
# variables or extra positional arguments into the job's own environment,
# so a conditional block here (`if [ "$SWEEP" = 1 ]`) silently fell through
# to the normal run below instead of triggering.

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

mkdir -p "$data"
mkdir -p "$data/$model"
mkdir -p "$data/$model/animation"

ll=0
numdata=1
for ((i=ll; i<ll+numdata; i++))
do
for m in 500.
do
   echo "calPzeta=$m"
   SIGMA=$m ./STOLAS $i
done
done
