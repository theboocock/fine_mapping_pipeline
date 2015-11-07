#!/usr/bin/env python
#
# Convert Z scores to pvalues
#

import argparse
from scipy.stats import norm


def print_p_values(zscores):
    """
        Read Zscores file, convert the Zscores to pvalues, and print them out
        along with the original data.
    """
    with open(zscores) as f:
        for line in f:
            split_line = line.split()
            rsid = line[0]
            zscore = line[1]
            print norm.pdf(zscore)
def main():
    parser = argparse.ArgumentParser(description="Make Zscores Pvalues")
    parser.add_argument("zscores", help="Zcores file")
    args = parser.parse_args()
    print_p_values(args.zscores)

if __name__=="__main__":
    main()
