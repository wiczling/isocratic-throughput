#!/bin/bash -l

#SBATCH --job-name=launcher
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH -p batch
#SBATCH --time=05:00:00
#SBATCH --mem=1gb

CHUNK=10
TOTAL=273

for ((start=1; start<=TOTAL; start+=CHUNK))
do
    end=$((start+CHUNK-1))

    if [ $end -gt $TOTAL ]; then
        end=$TOTAL
    fi

    echo "Submitting array ${start}-${end}"

    JOBID=$(sbatch --parsable --array=${start}-${end}%10 job_laplace.sh)

    echo "Submitted job ${JOBID}"

    echo "Waiting for completion..."

    while squeue -j ${JOBID} 2>/dev/null | grep -q ${JOBID}
    do
        sleep 30
    done

    echo "Chunk ${start}-${end} finished"

done

echo "All jobs completed"