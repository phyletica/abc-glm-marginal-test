#! /bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -j oe
#PBS -l jobflags=ADVRES:jro0014_lab.56281

if [ -n "$PBS_JOBNAME" ]
then
    source "${PBS_O_HOME}/.bash_profile"
    cd "$PBS_O_WORKDIR"
    module load gcc/5.3.0
fi

aln_num=33
seed=362847361
RANDOM=$seed

for run_num in 1 2 3
do
    out_path="run_stepping_stone_alignment_${aln_num}_run_${run_num}.out"
    rng_seed=$RANDOM
    ./run_stepping_stone.sh -a $aln_num -r $run_num -s $rng_seed 1>"$out_path" 2>&1
done
