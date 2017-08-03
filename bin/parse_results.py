#! /usr/bin/env python

import sys
import os
import re
import math
import glob

import project_util


number_pattern_str = r"[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?"
glm_density_pattern_str = r"marginal\s+density:\s+(?P<ml>" + number_pattern_str + ")"
glm_density_pattern = re.compile(glm_density_pattern_str, re.VERBOSE)


def get_max_abs_diff(values):
    mx = float("-inf")
    for i in range(len(values) - 1):
        for j in range(i + 1, len(values)):
            diff = math.fabs(values[i] - values[j])
            if diff > mx:
                mx = diff
    return mx

def parse_results(
        number_of_simulations = 100):
    if not os.path.exists(project_util.RESULTS_DIR):
        os.mkdir(project_util.RESULTS_DIR)
    results_path = os.path.join(project_util.RESULTS_DIR,
            "marginal-likelihood-results.txt")
    ss = []
    glm = []
    with open(results_path, 'w') as out:
        out.write("mean_ss_ml\tglm_ml\tmax_abs_diff_ss_ml\n")
        for i in range(100):
            abctoolbox_out_path = os.path.join(project_util.OBSERVED_DIR,
                    "sim-observed-alignment-{0}-abctoolbox.cfg.out".format(i+1))
            ss_ml_paths = glob.glob(os.path.join(project_util.OBSERVED_DIR,
                    "sim-observed-alignment-{0}-rb-ml-run-*.txt").format(i+1))
            if not len(ss_ml_paths) == 3:
                sys.stderr.write("found {0} rb ml paths for sim {1}\n".format(
                        len(ss_ml_paths), i + 1))
                continue
            assert len(ss_ml_paths) == 3, "found {0} rb ml paths for sim {1}".format(
                    len(ss_ml_paths), i + 1)
            ss_mls = []
            for path in ss_ml_paths:
                with open(path, 'r') as stream:
                    ss_mls.append(float(stream.read().strip()))
            mean_ss_ml = sum(ss_mls) / len(ss_mls)
            max_diff_ss_ml = get_max_abs_diff(ss_mls)
            with open(abctoolbox_out_path, 'r') as stream:
                glm_ml_matches = glm_density_pattern.findall(stream.read())
                assert len(glm_ml_matches) == 1
                glm_ml = float(glm_ml_matches[0])
            out.write("{mean_ss_ml}\t{glm_ml}\t{max_abs_diff_ss_ml}\n".format(
                    mean_ss_ml = mean_ss_ml,
                    glm_ml = glm_ml,
                    max_abs_diff_ss_ml = max_diff_ss_ml))
            ss.append(mean_ss_ml)
            glm.append(glm_ml)
    return ss, glm

def main_cli(argv = sys.argv):
    ss_mls, glm_mls = parse_results()


if __name__ == "__main__":
    main_cli()
