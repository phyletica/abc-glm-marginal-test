#! /usr/bin/env python

import os
import sys

BIN_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT_DIR = os.path.dirname(BIN_DIR)
OBSERVED_DIR = os.path.join(PROJECT_DIR, 'sim-observed-alignments')
PRIOR_DIR = os.path.join(PROJECT_DIR, 'abc-prior')

def main():
    sys.stdout.write("{0}\n".format(PROJECT_DIR))

if __name__ == '__main__':
    main()

