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
        out.write("true_edge_length\t"
                  "glm_ln_bf\t"
                  "trapezoidal_10000_ln_bf\t"
                  "correct_model_glm_ln_ml\t"
                  "correct_model_rectangular_100_ln_ml\t"
                  "correct_model_rectangular_1000_ln_ml\t"
                  "correct_model_rectangular_10000_ln_ml\t"
                  "correct_model_trapezoidal_100_ln_ml\t"
                  "correct_model_trapezoidal_1000_ln_ml\t"
                  "correct_model_trapezoidal_10000_ln_ml\t"
                  "correct_model_glm_mode_edge_length\t"
                  "correct_model_mcmc_mean_edge_length\t"
                  "vague_model_glm_ln_ml\t"
                  "vague_model_rectangular_100_ln_ml\t"
                  "vague_model_rectangular_1000_ln_ml\t"
                  "vague_model_rectangular_10000_ln_ml\t"
                  "vague_model_trapezoidal_100_ln_ml\t"
                  "vague_model_trapezoidal_1000_ln_ml\t"
                  "vague_model_trapezoidal_10000_ln_ml\t"
                  "vague_model_glm_mode_edge_length\t"
                  "vague_model_mcmc_mean_edge_length\n"
                  )
        correct_model_quadrature_diffs = []
        vague_model_quadrature_diffs = []
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
            trapezoidal_100_ln_ml = parse_value(trap_100_path)
            trapezoidal_1000_ln_ml = parse_value(trap_1000_path)
            trapezoidal_10000_ln_ml = parse_value(trap_10000_path)
            glm_mode_edge_length = parse_value(glm_mode_path)
            mcmc_mean_edge_length = parse_value(mcmc_mean_path)

            mx_diff = get_max_abs_diff((
                    # rectangular_100_ln_ml,
                    rectangular_1000_ln_ml,
                    rectangular_10000_ln_ml,
                    # trapezoidal_100_ln_ml,
                    trapezoidal_1000_ln_ml,
                    trapezoidal_10000_ln_ml))
            correct_model_quadrature_diffs.append(mx_diff)
            sys.stdout.write("\nmax diff for correct model: {0}\n".format(mx_diff))

            path_prefix += "vague-model-"
            vague_glm_ml_path     = path_prefix + "abctoolbox-glm-ml-estimate.txt"
            vague_rect_100_path   = path_prefix + "rectangular-ml-estimate-100.txt"
            vague_rect_1000_path  = path_prefix + "rectangular-ml-estimate-1000.txt"
            vague_rect_10000_path = path_prefix + "rectangular-ml-estimate-10000.txt"
            vague_trap_100_path   = path_prefix + "trapezoidal-ml-estimate-100.txt"
            vague_trap_1000_path  = path_prefix + "trapezoidal-ml-estimate-1000.txt"
            vague_trap_10000_path = path_prefix + "trapezoidal-ml-estimate-10000.txt"
            vague_glm_mode_path   = path_prefix + "abctoolbox-glm-edge-mode.txt"
            vague_mcmc_mean_path  = path_prefix + "mcmc-edge-mean.txt"

            vague_glm_ln_ml = math.log(parse_value(vague_glm_ml_path))
            vague_rectangular_100_ln_ml = parse_value(vague_rect_100_path)
            vague_rectangular_1000_ln_ml = parse_value(vague_rect_1000_path)
            vague_rectangular_10000_ln_ml = parse_value(vague_rect_10000_path)
            vague_trapezoidal_100_ln_ml = parse_value(vague_trap_100_path)
            vague_trapezoidal_1000_ln_ml = parse_value(vague_trap_1000_path)
            vague_trapezoidal_10000_ln_ml = parse_value(vague_trap_10000_path)
            vague_glm_mode_edge_length = parse_value(vague_glm_mode_path)
            vague_mcmc_mean_edge_length = parse_value(vague_mcmc_mean_path)

            mx_diff = get_max_abs_diff((
                    # vague_rectangular_100_ln_ml,
                    vague_rectangular_1000_ln_ml,
                    vague_rectangular_10000_ln_ml,
                    # vague_trapezoidal_100_ln_ml,
                    vague_trapezoidal_1000_ln_ml,
                    vague_trapezoidal_10000_ln_ml))
            vague_model_quadrature_diffs.append(mx_diff)
            sys.stdout.write("max diff for vague model: {0}\n".format(mx_diff))

            glm_ln_bf = glm_ln_ml - vague_glm_ln_ml
            trapezoidal_10000_ln_bf = trapezoidal_10000_ln_ml - vague_trapezoidal_10000_ln_ml
            out.write("{true_edge_length}\t"
                      "{glm_ln_bf}\t"
                      "{trapezoidal_10000_ln_bf}\t"
                      "{correct_model_glm_ln_ml}\t"
                      "{correct_model_rectangular_100_ln_ml}\t"
                      "{correct_model_rectangular_1000_ln_ml}\t"
                      "{correct_model_rectangular_10000_ln_ml}\t"
                      "{correct_model_trapezoidal_100_ln_ml}\t"
                      "{correct_model_trapezoidal_1000_ln_ml}\t"
                      "{correct_model_trapezoidal_10000_ln_ml}\t"
                      "{correct_model_glm_mode_edge_length}\t"
                      "{correct_model_mcmc_mean_edge_length}\t"
                      "{vague_model_glm_ln_ml}\t"
                      "{vague_model_rectangular_100_ln_ml}\t"
                      "{vague_model_rectangular_1000_ln_ml}\t"
                      "{vague_model_rectangular_10000_ln_ml}\t"
                      "{vague_model_trapezoidal_100_ln_ml}\t"
                      "{vague_model_trapezoidal_1000_ln_ml}\t"
                      "{vague_model_trapezoidal_10000_ln_ml}\t"
                      "{vague_model_glm_mode_edge_length}\t"
                      "{vague_model_mcmc_mean_edge_length}\n".format(
                          true_edge_length = true_edge_length,
                          glm_ln_bf = glm_ln_bf,
                          trapezoidal_10000_ln_bf = trapezoidal_10000_ln_bf,
                          correct_model_glm_ln_ml = glm_ln_ml,
                          correct_model_rectangular_100_ln_ml = rectangular_100_ln_ml,
                          correct_model_rectangular_1000_ln_ml = rectangular_1000_ln_ml,
                          correct_model_rectangular_10000_ln_ml = rectangular_10000_ln_ml,
                          correct_model_trapezoidal_100_ln_ml = trapezoidal_100_ln_ml,
                          correct_model_trapezoidal_1000_ln_ml = trapezoidal_1000_ln_ml,
                          correct_model_trapezoidal_10000_ln_ml = trapezoidal_10000_ln_ml,
                          correct_model_glm_mode_edge_length = glm_mode_edge_length,
                          correct_model_mcmc_mean_edge_length = mcmc_mean_edge_length,
                          vague_model_glm_ln_ml = vague_glm_ln_ml,
                          vague_model_rectangular_100_ln_ml = vague_rectangular_100_ln_ml,
                          vague_model_rectangular_1000_ln_ml = vague_rectangular_1000_ln_ml,
                          vague_model_rectangular_10000_ln_ml = vague_rectangular_10000_ln_ml,
                          vague_model_trapezoidal_100_ln_ml = vague_trapezoidal_100_ln_ml,
                          vague_model_trapezoidal_1000_ln_ml = vague_trapezoidal_1000_ln_ml,
                          vague_model_trapezoidal_10000_ln_ml = vague_trapezoidal_10000_ln_ml,
                          vague_model_glm_mode_edge_length = vague_glm_mode_edge_length,
                          vague_model_mcmc_mean_edge_length = vague_mcmc_mean_edge_length))

        correct_model_max_quadrature_diff = max(correct_model_quadrature_diffs)
        vague_model_max_quadrature_diff = max(vague_model_quadrature_diffs)
        sys.stdout.write("\nOverall max quadrature difference for correct model: {0}\n".format(
                correct_model_max_quadrature_diff))
        sys.stdout.write("Overall max quadrature difference for vague model: {0}\n".format(
                vague_model_max_quadrature_diff))

def main_cli(argv = sys.argv):
    parse_results()


if __name__ == "__main__":
    main_cli()
