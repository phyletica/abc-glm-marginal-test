#! /bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=40:00:00
#PBS -j oe
#PBS -l jobflags=ADVRES:jro0014_lab.56281

if [ -n "$PBS_JOBNAME" ]
then
    source "${PBS_O_HOME}/.bash_profile"
    cd "$PBS_O_WORKDIR"
fi

./calc_sum_stats.py 1>calc_sum_stats.sh.out 2>&1
