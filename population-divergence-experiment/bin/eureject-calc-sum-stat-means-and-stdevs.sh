#! /bin/sh

eureject -f ../sim-observed-alignments/sim-observed-alignment-1-sum-stats.txt \
        -k 0 \
        -n 100000 \
        -o ../abc-prior/sum-stat-means-and-stdevs.txt \
        ../abc-prior/prior-matrix.txt
