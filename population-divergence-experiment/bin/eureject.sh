#! /bin/bash

for observed_ss_file in ../sim-observed-alignments/sim-observed-alignment-*-sum-stats.txt
do
    posterior_path="${observed_ss_file/-sum-stats/-posterior-sample}"
    eureject -f "$observed_ss_file" -k 1000 -n 100000 -s ../abc-prior/sum-stat-means-and-stdevs.txt ../abc-prior/prior-matrix.txt > "$posterior_path"
done
