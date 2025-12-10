#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V
#$ -q all.q
#$ -N chaotic
#$ -o job_out
#$ -e job_out
#$ -pe OpenMP 14

#source /opt/intel/bin/compilervars.sh intel64
#export OMP_NUM_THREADS=$NSLOTS


# data=$(grep 'const std::string sdatadir' parameters.hpp | sed 's|//.*||' | sed 's/.*= *"\(.*\)".*/\1/')
# model=$(grep 'const std::string model' model.hpp | sed 's|//.*||' | sed 's/.*= *"\(.*\)".*/\1/')

# mkdir -p "$data"
# mkdir -p "$data""/$model"
# mkdir -p "$data""/$model/animation"

MODEL=$(grep '^#define MODEL' parameters.hpp | awk '{print $3}')

data=$(grep 'const std::string sdatadir' parameters.hpp \
        | sed 's|//.*||' \
        | sed 's/.*= *"\(.*\)".*/\1/')

if [ "$MODEL" -eq 0 ]; then
    model="chaotic"
elif [ "$MODEL" -eq 1 ]; then
    model="Starobinsky"
elif [ "$MODEL" -eq 2 ]; then
    model="USR"
else
    echo "Unknown MODEL=$MODEL"
    exit 1
fi

mkdir -p "$data"
mkdir -p "$data/$model"
mkdir -p "$data/$model/animation"

echo "/$model"

for ((i=0; i<1; i++))
do
    echo "Number $i"
	./STOLAS $i
done


