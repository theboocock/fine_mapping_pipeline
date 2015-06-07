# Copyright (c) 2015 Boocock James <james.boocock@otago.ac.nz>
# Author: Boocock James <james.boocock@otago.ac.nz>
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


__CAVIAR_BF_TEMPLATE__ = """
    caviarbf -z {0} -r {1} -t 0 -a {2} -n {3} -c {4} -o {5}
"""
__MODEL_SEARCH_TEMPLATE__ = """
    model_search -p 0 -m {0} -i {1] -o {2}
"""

import os
import logging

from fine_mapping_pipeline.expections.error_codes import *
from fine_mapping_pipeline.utils.shell import run_command

def run_caviarbf(output_directory, input_directory, sample_size, causal_snp_number=3, prior=0.1281429):
    """
        Runs caviar bf on the command-line     
    """
    logging.info("Starting to run caviarbf for each prepared input file")
    try:
        os.mkdir(output_directory)
    except OSError:
        logging.warning("Could not create output directory, it may exist already")
    input_files = os.path.join(input_directory, 'input.files')
    inputs = []
    with open(input_files) as in_f:
        for line in in_f:
            logging.info(line)
            inputs.append(line.strip())
    output_bfs = []  
    for input_f in inputs:
        z_score = os.path.join(input_directory,  input_f)
        ld = os.path.join(input_directory, input_f + '.Z')
        output_bf = os.path.join(output_directory, input_f +'.bf')
        output_bfs.append(output_bf)
        command = __CAVIAR_BF_TEMPLATE__.format(z_score, ld, prior, sample_size,
                                                causal_snp_number, output_bf)
        logging.debug("Running caviarbf with the command: {0}".format(output_bf))
        run_command(command, error=FAILED_CAVIARBF_RUN)
    logging.info("Running model_search on all bayes factor files")
    for input_f, output_bf in zip(input_files, output_bfs):
        no_snps = 0
        with open(os.path.join(output_directory, input_f)) as in_f:
            for no_snps, line in enumerate(in_f):
                pass
        output_basename = os.path.join(output_directory, input_f + '_caviar_run')
        command = __MODEL_SEARCH_TEMPLATE__.format(no_snps, input_f, output_basename)
        run_command(command, error=FAILED_CAVIARBF_RUN)
    logging.info("Caviar BF Run completed successfully")

def run_caviarbf_wrap(args):
    """
        Wraps the command-line arguments for running CaviarBF
    """
    _prepare_output_dir
