#! /usr/bin/env python

import os
import sys

BIN_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.dirname(BIN_DIR)
QSUB_DIR = os.path.join(BIN_DIR, "qsub-scripts")
MAIN_SCRIPT = os.path.join(BIN_DIR, "seq-divergence-model.py")
OUTPUT_DIR = os.path.join(PROJECT_DIR, "output")
RESULTS_DIR = os.path.join(PROJECT_DIR, "results")

def main():
    sys.stdout.write("{0}\n".format(PROJECT_DIR))

if __name__ == '__main__':
    main()

