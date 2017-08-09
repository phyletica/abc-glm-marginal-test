#! /bin/bash

current_dir="$(pwd)"

cd ../sim-observed-alignments

for abctoolbox_cfg in sim-observed-alignment-*-abctoolbox.cfg
do
    echo "$abctoolbox_cfg"
    ~/software/dev/PyMsBayes/bin/linux/ABCestimator "$abctoolbox_cfg" 1>"${abctoolbox_cfg}.out" 2>&1
done

cd "$current_dir"

