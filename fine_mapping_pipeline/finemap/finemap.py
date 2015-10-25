# Copyright (c) 2015 Boocock James <james.boocock@otago.ac.nz>
# Author: Boocock James <james.boocock@otago.ac.nz>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

import logging
import os
import tempfile
import sys
import numpy


from fine_mapping_pipeline.utils.shell import run_command
from fine_mapping_pipeline.expections.error_codes import *



__FINEMAP_TEMPLATE__="""
        finemap --sss --in-files {0} --n-ind {1} 
"""


def _read_matrix(matrix_file):
    """
        Read numpy matrix using the loadtxt function.
    """
    matrix = numpy.loadtxt(matrix_file, dtype='float')
    return matrix

def _write_matrix(matrix, output_matrix):
    """
        Writes a matrix to a file, used by remove surrogates
        to generate a new matrix for use with finemap.
    """
    numpy.savetxt(output_matrix, matrix, delimiter=' ', newline='\n')   

def _write_zscores(z_scores, output_zscores):
    with open(output_zscores, 'w') as z_score_out:
        for line in z_scores:
            z_score_out.write(line + "\n")

def remove_surrogates(matrix_file, z_score_file, surrogates_file=None):
    """
       Process the LD matrix to remove surrogates. 
    """
    matrix = _read_matrix(matrix_file)
    (nrow, ncol) = matrix.shape
    remove_snps =[]
    logging.info(matrix.shape)
    for i in range(nrow):
        for j in range(i, ncol):
            if i != j:
                if matrix[i][j] == 1:
                    remove_snps.append(j)
    matrix = numpy.delete(matrix, remove_snps, axis=1)
    matrix = numpy.delete(matrix, remove_snps, axis=0)
    z_scores = [] 
    surrogates_output = surrogates_file
    if surrogates_output is not None:
        print surrogates_file
        surrogates_output = open(surrogates_file,'w')
    with open(z_score_file) as z:
        for i, line in enumerate(z):
            if i not in remove_snps:
                z_scores.append(line.strip())
            elif surrogates_output is not None:
                surrogates_output.write(line.split()[0] + "\n")
    return (matrix, z_scores)

def _write_k_file(output_k, causal_snp_number):
    """
        Writes the Casual SNP K file. 

        Let's just write the file that I used in the example.

        0.6 0.3 0.1
    """

    ## TODO Change this horrible default
    with open(output_k, 'w') as out:
        out.write("0.6 0.3 0.1\n")



def run_finemap(input_directory, causal_snp_number, output_directory, sample_size):
    """
       Runs the finemapping pipeline 

       TODO: Remove causal SNP number argument - not needed.
    """
    # Get the matrix
    try:
        os.mkdir(output_directory)
    except OSError:
        pass
    input_files = os.path.join(input_directory, 'input.files')
    inputs = []
    with open(input_files) as in_f:
        for line in in_f:
            logging.info(line)
            prefix = line.strip()
            matrix_file = os.path.join(input_directory, prefix +".matrix")
            z_score_file = os.path.join(input_directory, prefix + ".Z")
            output_surrogates = os.path.join(input_directory, prefix + ".surrogates")
            (matrix, z_scores) = remove_surrogates(matrix_file, z_score_file, output_surrogates)
            try:
                tmp_directory = tempfile.mkdtemp() 
            except OSError:
                logging.error('Could not create a temporary directory for likelihood testing')
                sys.exit(OS_ERROR)
            output_matrix = os.path.join(tmp_directory, prefix + ".ld")
            output_z = os.path.join(tmp_directory, prefix +".z")
            output_k = os.path.join(tmp_directory, prefix +".k")
            _write_matrix(matrix, output_matrix)
            _write_zscores(z_scores, output_z)
            _write_k_file(output_k, causal_snp_number)
            # Let's run the analysis
            finemap_cl = __FINEMAP_TEMPLATE__.format(os.path.join(tmp_directory,prefix), sample_size) 
            logging.info(finemap_cl)
            run_command(finemap_cl, error=FAILED_FINEMAP_RUN) 
    # TODO extract the output from the analysis and interpret the results :O

def run_finemap_wrap(args):
    """
        Wraps the run finemap function so that it can be used from the command line
    """
    input_directory = args.input_directory
    causal_snp_number = args.causal_snp_number
    number_of_ind = args.number_of_individuals
    output_directory = args.output_directory
    run_finemap(input_directory, causal_snp_number, output_directory, number_of_ind)



