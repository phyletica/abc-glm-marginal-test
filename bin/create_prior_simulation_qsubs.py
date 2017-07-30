#! /usr/bin/env python

import sys
import os
import random

import project_util

def get_pbs_header(restrict_nodes = False, walltime = "1:00:00"):
    s = ("#! /bin/sh\n"
         "#PBS -l nodes=1:ppn=1\n"
         "#PBS -l walltime={0}\n"
         "#PBS -j oe\n".format(walltime))
    if restrict_nodes:
        s += "#PBS -l jobflags=ADVRES:jro0014_lab.56281\n"
    s += ("\n"
          "if [ -n \"$PBS_JOBNAME\" ]\n"
          "then\n"
          "    source \"${PBS_O_HOME}/.bash_profile\"\n"
          "    cd \"$PBS_O_WORKDIR\"\n"
          "    module load gcc/5.3.0\n"
          "fi\n\n")
    return s

def main():
    pbs_header = get_pbs_header()
    rng = random.Random()
    rng.seed(951640981)
    number_of_batches = 100
    batch_size = 1000
    batch_script_path = os.path.join(project_util.BIN_DIR,
            "simulate_prior_alignments.sh")
    qsub_dir = os.path.join(project_util.BIN_DIR,
            "simulate_prior_alignments_batches")
    rel_batch_script_path = os.path.relpath(batch_script_path, qsub_dir)
    if not os.path.exists(qsub_dir):
        os.mkdir(qsub_dir)
    for i in range(number_of_batches):
        qsub_script_path = os.path.join(qsub_dir,
                "simulate_prior_batch_{0}.sh".format(i +  1))
        stdout_path = os.path.basename(qsub_script_path) + ".out"
        output_dir = os.path.join(project_util.PRIOR_DIR,
                "batch-{0}".format(i + 1))
        rel_output_dir = os.path.relpath(output_dir, qsub_dir)
        cmd_line = "{p} -n {n} -o {o} -s {s} 1>{stdout} 2>&1\n".format(
                p = rel_batch_script_path,
                n = batch_size,
                o = rel_output_dir,
                s = rng.randint(1, 999999999999),
                stdout = stdout_path)
        with open(qsub_script_path, "w") as out:
            out.write("{0}{1}".format(pbs_header, cmd_line))

if __name__ == "__main__":
    main()
