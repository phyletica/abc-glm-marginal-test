#! /usr/bin/env python

import sys
import os
import math
import glob

import project_util


def parse_value(path, header = False):
    with open(path, "r") as stream:
        if header:
            stream.next()
        value_str = stream.next().strip()
    return float(value_str)

def get_max_abs_diff(values):
    mx = float("-inf")
    for i in range(len(values) - 1):
        for j in range(i + 1, len(values)):
            diff = math.fabs(values[i] - values[j])
            if diff > mx:
                mx = diff
    return mx

def parse_results():
    if not os.path.exists(project_util.RESULTS_DIR):
        os.mkdir(project_util.RESULTS_DIR)
    results_path = os.path.join(project_util.RESULTS_DIR,
            "results.txt")
    with open(results_path, 'w') as out:
        out.write("glm_ln_ml\t"
                  "rectangular_100_ln_ml\t"
                  "rectangular_1000_ln_ml\t"
                  "rectangular_10000_ln_ml\t"
                  "trapezoidal_100_ln_ml\t"
                  "trapezoidal_1000_ln_ml\t"
                  "trapezoidal_10000_ln_ml\t"
                  "true_edge_length\t"
                  "glm_mode_edge_length\t"
                  "mcmc_mean_edge_length\n")
        for i in range(1, 101):
            path_prefix = os.path.join(project_util.OUTPUT_DIR, "sim-{0}-".format(i))
            true_path       = path_prefix + "true-parameter-value.txt"
            glm_ml_path     = path_prefix + "abctoolbox-glm-ml-estimate.txt"
            rect_100_path   = path_prefix + "rectangular-ml-estimate-100.txt"
            rect_1000_path  = path_prefix + "rectangular-ml-estimate-1000.txt"
            rect_10000_path = path_prefix + "rectangular-ml-estimate-10000.txt"
            trap_100_path   = path_prefix + "trapezoidal-ml-estimate-100.txt"
            trap_1000_path  = path_prefix + "trapezoidal-ml-estimate-1000.txt"
            trap_10000_path = path_prefix + "trapezoidal-ml-estimate-10000.txt"
            glm_mode_path   = path_prefix + "abctoolbox-glm-edge-mode.txt"
            mcmc_mean_path  = path_prefix + "mcmc-edge-mean.txt"

            true_edge_length = parse_value(true_path, header = True)
            glm_ln_ml = math.log(parse_value(glm_ml_path))
            rectangular_100_ln_ml = parse_value(rect_100_path)
            rectangular_1000_ln_ml = parse_value(rect_1000_path)
            rectangular_10000_ln_ml = parse_value(rect_10000_path)
            trapezoidal_100_ln_ml = parse_value(rect_100_path)
            trapezoidal_1000_ln_ml = parse_value(rect_1000_path)
            trapezoidal_10000_ln_ml = parse_value(rect_10000_path)
            glm_mode_edge_length = parse_value(glm_mode_path)
            mcmc_mean_edge_length = parse_value(mcmc_mean_path)

            mx_diff = get_max_abs_diff((
                    rectangular_100_ln_ml,
                    rectangular_1000_ln_ml,
                    rectangular_10000_ln_ml,
                    trapezoidal_100_ln_ml,
                    trapezoidal_1000_ln_ml,
                    trapezoidal_10000_ln_ml))
            sys.stdout.write("max diff: {0}\n".format(mx_diff))

            out.write("{glm_ln_ml}\t"
                      "{rectangular_100_ln_ml}\t"
                      "{rectangular_1000_ln_ml}\t"
                      "{rectangular_10000_ln_ml}\t"
                      "{trapezoidal_100_ln_ml}\t"
                      "{trapezoidal_1000_ln_ml}\t"
                      "{trapezoidal_10000_ln_ml}\t"
                      "{true_edge_length}\t"
                      "{glm_mode_edge_length}\t"
                      "{mcmc_mean_edge_length}\n".format(
                          glm_ln_ml = glm_ln_ml,
                          rectangular_100_ln_ml = rectangular_100_ln_ml,
                          rectangular_1000_ln_ml = rectangular_1000_ln_ml,
                          rectangular_10000_ln_ml = rectangular_10000_ln_ml,
                          trapezoidal_100_ln_ml = trapezoidal_100_ln_ml,
                          trapezoidal_1000_ln_ml = trapezoidal_1000_ln_ml,
                          trapezoidal_10000_ln_ml = trapezoidal_10000_ln_ml,
                          true_edge_length = true_edge_length,
                          glm_mode_edge_length = glm_mode_edge_length,
                          mcmc_mean_edge_length = mcmc_mean_edge_length))

def main_cli(argv = sys.argv):
    parse_results()


if __name__ == "__main__":
    main_cli()
