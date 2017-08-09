#! /bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=0:10:00
#PBS -j oe

if [ -n "$PBS_JOBNAME" ]
then
    source "${PBS_O_HOME}/.bash_profile"
    cd "$PBS_O_WORKDIR"
    module load gcc/5.3.0
fi

./simulate_prior_alignments.sh -n 1000 -o ../abc-prior/batch-16 -s 286229179 1>simulate_prior_alignments_batch_16.sh.out 2>&1
