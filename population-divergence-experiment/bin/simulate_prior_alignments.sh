#! /bin/bash

output_dir="."
number_of_sims=10
rng_seed=$RANDOM

usage () {
    echo ""
    echo "usage: $0 [-h|--help] [-n|--number-of-sims <N>] [-o|--output-dir <PATH>] [-s|--seed <N>]"
    echo "  -h|--help           Show help message and exit."
    echo "  -n|--number-of-sims Number of alignments to simulate."
    echo "                      Default: ${number_of_sims}."
    echo "  -o|--output-dir     Directory to write simulations."
    echo "                      Default: \"${output_dir}\"."
    echo "  -s|--seed           Seed for RevBayes random number generator."
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
        -n| --number-of-sims)
            shift
            number_of_sims="$1"
            ;;
        -o| --output-dir)
            shift
            output_dir="$1"
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
      output_dir = \"${output_dir}\";
      number_of_sims = ${number_of_sims};
      source(\"simulate_prior_alignments.rev\");"

echo $args | rb
