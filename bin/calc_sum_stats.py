#! /usr/bin/env python

import sys
import os
import glob

import seqsift
from seqsift.seqops import seqstats

import project_util


def calc_sum_stats(fasta_path, pop1_prefix = "sp1", pop2_prefix = "sp2"):
    seq_iter = seqsift.dataio.get_seq_iter_from_file(
            fasta_path,
            format = "fasta")
    pop1_seqs = []
    pop2_seqs = []
    nsites = None
    for seq in seq_iter:
        if nsites is None:
            nsites = len(seq.seq)
        else:
            assert len(seq.seq) == nsites
        if seq.id.startswith(pop1_prefix):
            pop1_seqs.append(seq)
        elif seq.id.startswith(pop2_prefix):
            pop2_seqs.append(seq)
        else:
            raise Exception("Unexpected taxon label {0} in {1}".format(
                seq.id, fasta_path))
    assert len(pop1_seqs) == len(pop2_seqs)

    return seqstats.get_population_pair_diversity_summary(
            pop1_seqs,
            pop2_seqs,
            per_site = True,
            aligned = True)

def calc_sum_stats_for_observed_alignments(number_of_alignments = 10):
    for i in range(number_of_alignments):
        observed_path = os.path.join(project_util.OBSERVED_DIR,
                "sim-observed-alignment-{0}.fasta".format(i + 1))
        path_prefix = os.path.splitext(observed_path)[0]
        sum_stat_path = path_prefix + "-sum-stats.txt"
        sum_stats = calc_sum_stats(observed_path)
        with open(sum_stat_path, "w") as out:
            out.write("pi_net\tpi_within\n")
            out.write("{pi_net}\t{pi_within}\n".format(**sum_stats))

def calc_sum_stats_for_prior_alignments():
    batch_dirs = sorted(glob.glob(os.path.join(
            project_util.PRIOR_DIR,
            "batch*")))
    prior_matrix_path = os.path.join(project_util.PRIOR_DIR,
            "prior-matrix.txt")
    with open(prior_matrix_path, "w") as out:
        out.write("root_height\tpop_size\tpi_net\tpi_within\n")
        for batch_dir in batch_dirs:
            parameter_values_path = os.path.join(batch_dir,
                    "parameter-values.txt")
            with open(parameter_values_path, "r") as stream:
                header = stream.next().strip()
                assert header == "root_height\tpop_size"
                for i, line in enumerate(stream):
                    l = line.strip()
                    alignment_path = os.path.join(batch_dir,
                            "prior-sample-alignment-{0}.fasta".format(i + 1))
                    sum_stats = calc_sum_stats(alignment_path)
                    l += "\t{pi_net}\t{pi_within}".format(**sum_stats)
                    out.write("{0}\n".format(l))

def main_cli(argv = sys.argv):
    seqsift_info_path = os.path.join(project_util.BIN_DIR,
            "seqsift-version-info.txt")
    with open(seqsift_info_path, "w") as out:
        out.write("{0}\n".format(seqsift.get_description()))
    calc_sum_stats_for_observed_alignments(10)
    calc_sum_stats_for_prior_alignments()


if __name__ == "__main__":
    main_cli()
