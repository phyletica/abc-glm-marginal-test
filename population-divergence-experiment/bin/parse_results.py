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

def parse_revbayes_ml_estimates(ss_ml_paths):
    ss_mls = []
    for path in ss_ml_paths:
        with open(path, 'r') as stream:
            ss_mls.append(float(stream.read().strip()))
    return ss_mls

def parse_posterior_sample(path):
    heights = []
    sizes = []
    with open(path, "r") as stream:
        header = stream.next().strip().split()
        assert header[0:2] == ["root_height", "pop_size"]
        for line in stream:
            h, s = (float(x) for x in line.strip().split()[0:2])
            heights.append(h)
            sizes.append(s)
    assert len(heights) == len(sizes)
    mean_height = sum(heights) / len(heights)
    mean_size = sum(sizes) / len(sizes)
    return mean_height, mean_size

def parse_glm_estimates(path):
    with open(path, "r") as stream:
        header = stream.next().strip().split()
        assert header == ["what", "root_height", "pop_size"]
        modes = stream.next().strip().split()
        assert modes[0] == "mode"
        assert len(modes) == 3
    return float(modes[1]), float(modes[2])

def parse_abctoolbox_stdout(path):
    with open(path, "r") as stream:
        glm_ml_matches = glm_density_pattern.findall(stream.read())
        assert len(glm_ml_matches) == 1
    return float(glm_ml_matches[0])

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
            "marginal-likelihood-results.txt")
    true_parameter_values_path = os.path.join(project_util.OBSERVED_DIR,
            "parameter-values.txt")
    with open(results_path, 'w') as out:
        out.write("mean_ss_ml\t"
                  "glm_ml\t"
                  "max_abs_diff_ss_ml\t"
                  "true_height\t"
                  "post_sample_mean_height\t"
                  "glm_mode_height\t"
                  "true_size\t"
                  "post_sample_mean_size\t"
                  "glm_mode_size\n")
        with open(true_parameter_values_path, "r") as true_parameter_stream:
            true_parameter_header = true_parameter_stream.next().strip().split()
            assert true_parameter_header == ["root_height", "pop_size"]
            for i, true_parameter_line in enumerate(true_parameter_stream):
                true_height, true_size = (float(x) for x in true_parameter_line.strip().split())

                posterior_sample_path = os.path.join(project_util.OBSERVED_DIR,
                        "sim-observed-alignment-{0}-posterior-sample.txt".format(i+1))
                abctoolbox_out_path = os.path.join(project_util.OBSERVED_DIR,
                        "sim-observed-alignment-{0}-abctoolbox.cfg.out".format(i+1))
                glm_parameter_estimates_path = os.path.join(project_util.OBSERVED_DIR,
                        "sim-observed-alignment-{0}-abc-glm-PosteriorCharacteristics_Obs0.txt".format(i+1))
                ss_ml_paths = glob.glob(os.path.join(project_util.OBSERVED_DIR,
                        "sim-observed-alignment-{0}-rb-ml-run-*.txt").format(i+1))
                assert len(ss_ml_paths) == 3, "found {0} rb ml paths for sim {1}".format(
                        len(ss_ml_paths), i + 1)

                ss_mls = parse_revbayes_ml_estimates(ss_ml_paths)
                mean_ss_ml = sum(ss_mls) / len(ss_mls)
                max_diff_ss_ml = get_max_abs_diff(ss_mls)

                glm_ml = parse_abctoolbox_stdout(abctoolbox_out_path) 
                post_mean_height, post_mean_size = parse_posterior_sample(
                        path = posterior_sample_path)
                glm_mode_height, glm_mode_size = parse_glm_estimates(
                        path = glm_parameter_estimates_path)
                out.write("{mean_ss_ml}\t"
                          "{glm_ml}\t"
                          "{max_abs_diff_ss_ml}\t"
                          "{true_height}\t"
                          "{post_sample_mean_height}\t"
                          "{glm_mode_height}\t"
                          "{true_size}\t"
                          "{post_sample_mean_size}\t"
                          "{glm_mode_size}\n".format(
                        mean_ss_ml = mean_ss_ml,
                        glm_ml = glm_ml,
                        max_abs_diff_ss_ml = max_diff_ss_ml,
                        true_height = true_height,
                        post_sample_mean_height = post_mean_height,
                        glm_mode_height = glm_mode_height,
                        true_size = true_size,
                        post_sample_mean_size = post_mean_size,
                        glm_mode_size = glm_mode_size
                        ))

def main_cli(argv = sys.argv):
    parse_results()


if __name__ == "__main__":
    main_cli()
