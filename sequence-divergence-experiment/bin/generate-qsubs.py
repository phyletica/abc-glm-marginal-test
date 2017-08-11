#! /usr/bin/env python

import sys
import os
import random

import project_util

_RNG = random.Random()

def get_pbs_header(restrict_nodes = True, walltime = "5:00:00"):
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
          "    which python\n"
          "fi\n\n")
    return s

def main():
    restrict_nodes = True
    walltime = "5:00:00"

    seq_length = 10000
    abc_prior_samples = 100000
    abc_posterior_samples = 5000
    mcmc_generations = 10000
    mcmc_sample_frequency = 10
    prior_lower = 1.0 / seq_length
    prior_upper = 0.1
    vague_prior_lower = 1.0 / (10 * seq_length)
    vague_prior_upper = 0.15

    if not os.path.exists(project_util.QSUB_DIR):
        os.mkdir(project_util.QSUB_DIR)
    if not os.path.exists(project_util.OUTPUT_DIR):
        os.mkdir(project_util.OUTPUT_DIR)
    relative_output_dir = os.path.relpath(project_util.OUTPUT_DIR,
            project_util.QSUB_DIR)
    relative_script_path = os.path.relpath(project_util.MAIN_SCRIPT,
            project_util.QSUB_DIR)

    rng = random.Random()
    rng.seed(1234)

    number_of_simulations = 100
    for i in range(number_of_simulations):
        qsub_path = os.path.join(project_util.QSUB_DIR,
                "sim-{0}.sh".format(i + 1))
        stdout_path = os.path.basename(qsub_path) + ".out"

        output_prefix = os.path.join(relative_output_dir,
                "sim-{0}".format(i + 1))
        seed = rng.randint(1, 999999999)

        with open(qsub_path, "w") as out:
            out.write(get_pbs_header(
                    restrict_nodes = restrict_nodes,
                    walltime = walltime))
            cmd = ("python {relative_script_path} --seq-length {seq_length} "
                   "--abc-prior-samples {abc_prior_samples} "
                   "--abc-posterior-samples {abc_posterior_samples} "
                   "--mcmc-generations {mcmc_generations} "
                   "--mcmc-sample-frequency {mcmc_sample_frequency} "
                   "--prior-lower {prior_lower} "
                   "--prior-upper {prior_upper} "
                   "--vague-prior-lower {vague_prior_lower} "
                   "--vague-prior-upper {vague_prior_upper} "
                   "--output-prefix {output_prefix} "
                   "{seed} 1>{stdout_path} 2>&1\n".format(
                            relative_script_path = relative_script_path,
                            seq_length = seq_length,
                            abc_prior_samples = abc_prior_samples,
                            abc_posterior_samples = abc_posterior_samples,
                            mcmc_generations = mcmc_generations,
                            mcmc_sample_frequency = mcmc_sample_frequency,
                            prior_lower = prior_lower,
                            prior_upper = prior_upper,
                            vague_prior_lower = vague_prior_lower,
                            vague_prior_upper = vague_prior_upper,
                            output_prefix = output_prefix,
                            seed = seed,
                            stdout_path = stdout_path))
            out.write(cmd)

if __name__ == "__main__":
    main()
