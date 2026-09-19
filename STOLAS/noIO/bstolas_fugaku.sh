#!/bin/bash

#PJM -L "node=1"
#PJM -L "rscgrp=int"
#PJM -L "elapse=05:00:00"
#PJM -g hp260514
#PJM -x PJM_LLIO_GFSCACHE=/vol0005
#PJM -j
#PJM -S
#PJM -o job_out

module load lang/tcsds-1.2.43

export OMP_NUM_THREADS=48

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
