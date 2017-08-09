#! /usr/bin/env python

import sys
import os
import glob

import project_util

def get_abctoolbox_config_str(
        posterior_path,
        sum_stats_path,
        output_prefix,
        dirac_peak_width = 0.002,
        posterior_density_points = 500):
    cfg = ("estimationType standard\n"
           "simName {posterior_path}\n"
           "obsName {sum_stats_path}\n"
           "params 1-2\n"
           "diracPeakWidth {dirac_peak_width}\n"
           "posteriorDensityPoints {posterior_density_points}\n"
           "stadardizeStats 1\n"
           "maxReadSims 100000\n"
           "outputPrefix {output_prefix}\n".format(
                    posterior_path = posterior_path,
                    sum_stats_path = sum_stats_path,
                    output_prefix = output_prefix,
                    dirac_peak_width = dirac_peak_width,
                    posterior_density_points = posterior_density_points))
    return cfg

def main_cli(argv = sys.argv):
    posterior_paths = sorted(glob.glob(os.path.join(project_util.OBSERVED_DIR,
            "sim-observed-alignment-*-posterior-sample.txt")))

    for post_path in posterior_paths:
        post_file_name = os.path.basename(post_path)
        sum_stats_path = post_file_name.replace("posterior-sample", "sum-stats")
        output_prefix = post_file_name.replace("posterior-sample.txt", "abc-glm-")
        cfg_path = post_path.replace("posterior-sample.txt", "abctoolbox.cfg")
        cfg_str = get_abctoolbox_config_str(
                posterior_path = post_file_name,
                sum_stats_path = sum_stats_path,
                output_prefix = output_prefix,
                dirac_peak_width = 0.002,
                posterior_density_points = 500)
        with open(cfg_path, "w") as out:
            out.write(cfg_str)

if __name__ == "__main__":
    main_cli()
