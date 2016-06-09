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

from shutil import rmtree

from fine_mapping_pipeline.utils.shell import run_command
from fine_mapping_pipeline.expections.error_codes import *
from scipy.stats import chi2

__PAINTOR_TEMPLATE__="""
     /Users/smilefreak/PAINTOR_V3.0/PAINTOR -input {0}/input.files  -in {0}/ -out {1}/ -Zhead {2} -LDname {3} -c {4} -only-enrichment 
"""

def _do_lrt(null_lrt, annot_lrt):
    """
        Perform the likelihood ratio test using standard LRT. 
    """
    test = -2*(null_lrt - annot_lrt)
    test_result = 1 - chi2.cdf(test, 1)
    return test_result

def _get_likelihood(input_directory, i, annotation, causal_snp_number, null_likelihood, populations):
    """
        Use the likelhood ratio test to determine the sigfinance of an individual annotation for PAINTOR. 
    """
    try:
        output_directory = tempfile.mkdtemp() 
    except OSError:
        logging.error('Could not create a temporary directory for likelihood testing')
        sys.exit(OS_ERROR)
    logging.info("Testing annotation {0} for significance using the LRT".format(annotation))
    ld_name = populations.split(',')
    ld_suffix = ['LD.'+o for o in ld_name]
    ld_suffix = ",".join(ld_suffix)
    command = __PAINTOR_TEMPLATE__.format(input_directory, output_directory,
                                         populations, ld_suffix, causal_snp_number) 
    if i != -1:
        command += '-annotations ' + i 
    logging.info("Paintor command = {0}".format(command))
    run_command(command, error=FAILED_PAINTOR_RUN, exit_on_failure=False)
    likelihood = 0.0
    try:
        with open(os.path.join(output_directory,"Log.Likelihood")) as likeli_file:
            try:
                likelihood = float(likeli_file.readline())
            except ValueError:
                logging.error("Could not convert PAINTOR likelihood output to a float")
                sys.exit(FAILED_PAINTOR_RUN)
    except IOError:
        likelihood = null_likelihood
        logging.warning("No likelihood was produced for this annotation, paintor run must have failed")
    try:
        rmtree(output_directory)
    except OSError:
        logging.error("Could not remove a temporary directory created for the paintor Run")
        sys.exit(OS_ERROR)
    return likelihood

def _select_annotations(input_directory, annotation_header, 
                        causal_snp_number, populations, p_value_threshold=0.05):
    p_value_threshold = p_value_threshold
    best_annotations = []
    logging.info("Selecting annotations to use with the LRT")
    null_likelihood = _get_likelihood(input_directory, -1, "", causal_snp_number,0, populations)
    logging.debug("Null model likelihood: {0}".format(null_likelihood))
    while annotation_header:
        lrts = []
        for i, annotation in enumerate(annotation_header):
            logging.info("Testing annotation {0}".format(annotation))
            temp_likelhood = _get_likelihood(input_directory, annotation,
                                           annotation_header, causal_snp_number, null_likelihood, populations)
            logging.debug(temp_likelhood)
            lrt = _do_lrt(null_likelihood, temp_likelhood)
            lrts.append((i, lrts))
            logging.debug(lrt)
            logging.info("DETECTED: annotation {0}: pvalue = {1}".format(annotation, lrt))
    return best_annotations


def run_paintor(input_directory, output_directory, populations,
                auto_select_annotations=False, causal_snp_number=3):
    """
        Function runs PAINTOR and selections the annotations for using Downstream.
    """
    # TODO first thing is to test to see which regions fail using just a -c 1
    logging.info("Running Paintor")
    try:
        temp_output_directory = tempfile.mkdtemp() 
    except OSError:
        logging.error('Could not create a temporary directory for testing for PAINTOR failure')
        sys.exit(OS_ERROR)
    keep_these_regions = []
    loci = []
    with open(os.path.join(input_directory, 'input.files')) as f:
        for i, line in enumerate(f):
            if i ==0:
                with open(os.path.join(input_directory, line.strip() +".annotations")) as annotations_all:
                    header = annotations_all.readline()
                    annotation_header = header.strip().split()
            loci.append(line.strip())
    if auto_select_annotations:
        best_annotations = _select_annotations(input_directory, annotation_header, causal_snp_number, populations)
    else:
        best_annotations = range(len(header_line))
    command = __PAINTOR_TEMPLATE__.format(input_directory, output_directory,
                                          os.path.join(input_directory, 'input.files'),
                                          causal_snp_number)
    if len(best_annotations) > 0:
        command += '-i '+ ','.join([str(o) for o in best_annotations])
    run_command(command)
    logging.info("Header Annotation")
    logging.info([str(o) for i, o in enumerate(header_line) if i in best_annotations])
    logging.info("Finished running Paintor")

def run_paintor_wrap(args):
    """
        Wraps the run paintor function so that it can be used from the command line
    """

    auto_select_annotations = args.auto_select_annotations
    input_directory = args.input_directory
    output_directory = args.output_directory
    causal_snp_number = args.causal_snp_number
    populations = args.populations
    try:
        os.mkdir(output_directory)
    except OSError:
        pass
    run_paintor(input_directory, output_directory, populations, auto_select_annotations, causal_snp_number)
