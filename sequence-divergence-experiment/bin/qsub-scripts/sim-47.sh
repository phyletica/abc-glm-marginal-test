#! /bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=5:00:00
#PBS -j oe
#PBS -l jobflags=ADVRES:jro0014_lab.56281

if [ -n "$PBS_JOBNAME" ]
then
    source "${PBS_O_HOME}/.bash_profile"
    cd "$PBS_O_WORKDIR"
    which python
fi

python ../seq-divergence-model.py --seq-length 10000 --abc-prior-samples 100000 --abc-posterior-samples 5000 --mcmc-generations 10000 --mcmc-sample-frequency 10 --prior-lower 0.0001 --prior-upper 0.1 --vague-prior-lower 1e-05 --vague-prior-upper 0.15 --output-prefix ../../output/sim-47 144245722 1>sim-47.sh.out 2>&1
