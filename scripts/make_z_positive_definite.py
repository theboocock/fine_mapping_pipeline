#!/usr/bin/env python
#
# Take a correlation matrix that has been used for paintor analysis and convert into a positive definite matrix
#


import argparse
import numpy

def main():
    """
        Take a correlation matrix into python and convert the matrix into positive definite.
    """

    parser = argparse.ArgumentParser(description="Correlation matrix to positive definite")
    parser.add_argument("correlation_matrix")
    parser.add_argument("-n","--no-snps",dest="number_of_snps")
    args = parser.parse_args()
    # Todo replace results here.
    n = int(args.number_of_snps)
    cor_mat = numpy.loadtxt(args.correlation_matrix, dtype=numpy.float32)
    cor_mat =  numpy.power(cor_mat, 2)
    print(cor_mat)
    ## Write out matrix

if __name__ == "__main__":
    main()

