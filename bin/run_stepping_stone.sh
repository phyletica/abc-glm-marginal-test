#! /bin/bash

aln_num=1
run_num=1
rng_seed=$RANDOM

usage () {
    echo ""
    echo "usage: $0 [-h|--help] [-n|--number-of-sims <N>] [-o|--output-dir <PATH>] [-s|--seed <N>]"
    echo "  -h|--help       Show help message and exit."
    echo "  -a|--aln-num    Number for observed alignment to analyze."
    echo "                      Default: ${aln_num}."
    echo "  -r|--run-num    Run number."
    echo "                      Default: ${run_num}."
    echo "  -s|--seed       Seed for RevBayes random number generator."
    echo "                      Default: Use \$RANDOM."
    echo ""
}

# process args
if [ "$(echo "$@" | grep -c "=")" -gt 0 ]
then
    echo "ERROR: Do not use '=' for arugments. For example, use"
    echo "'--nthreads 2' instead of '--nthreads=2'."
    exit 1
fi

extra_args=""
while [ "$1" != "" ]
do
    case $1 in
        -h| --help)
            usage
            exit
            ;;
        -a| --aln-num)
            shift
            aln_num="$1"
            ;;
        -r| --run-num)
            shift
            run_num="$1"
            ;;
        -s| --seed)
            shift
            rng_seed="$1"
            ;;
        * )
            extra_args="$extra_args $1"
    esac
    shift
done

if [ -n "$extra_args" ]
then
    echo "ERROR: unexpected arguments: $extra_args"
    usage
    exit 1
fi

args="rng_seed = ${rng_seed};
      aln_num = ${aln_num};
      run_num = ${run_num};
      source(\"run_stepping_stone.rev\");"

echo $args | rb

rm "../sim-observed-alignments/sim-observed-alignment-${aln_num}-rb-pp-run-${run_num}_stone_*.log"
