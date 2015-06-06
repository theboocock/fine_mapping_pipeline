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
import sys

from fine_mapping_pipeline.expections.error_codes import *
from fine_mapping_pipeline.utils.shell import *

__VCF_TO_PLINK_TEMPLATE__ = """
    plink --vcf {0} --recode --out {1} 
"""
__PLINK_TO_LD_MATRIX__ = """
    plink --file {0} --matrix --out {1} --r --allow-no-sex
"""

def _add_dimensions_to_file(locus_f):
    """
        Function adds the dimensions for the LD file so it can be used in Paintor
    """
    ld_lines = []
    i = 0
    with open(locus_f) as ld_file:
        for i, line in enumerate(ld_file):
            ld_lines.append(line)
    no_lines = i + 1
    file_out = locus_f.split('.matrix')[0] + '.LD'
    with open(file_out, 'w' ) as paintor_ld:
        paintor_ld.write(str(no_lines) + ' ' + str(no_lines) + '\n')
        for line in ld_lines:
            paintor_ld.write(line)

def vcf_to_plink(locus, output_directory ,vcf):
    """
        Uses PLINK to convert a 1000 genomes VCF file to plink format 
    """
    logging.info("Converting {0} to PLINK format".format(vcf))
    command = __VCF_TO_PLINK_TEMPLATE__.format(vcf, locus)
    run_command(command)
    try:
        #Rename
        os.rename(locus +'.map', os.path.join(output_directory, locus + '.map'))
        os.rename(locus +'.ped', os.path.join(output_directory, locus + '.ped'))
        #Remove
        os.remove(locus + '.log')
        os.remove(locus + '.nosex')
    except OSError as e:
        logging.error("Could not move PLINK files {0}".format(e))
        sys.exit(OS_ERROR)

def _remove_plink_files(output_directory, locus):
    """
        Remove the plinkfile after creating the basename file.
    """
    ped = os.path.join(output_directory, locus + '.ped')
    map_f = os.path.join(output_directory, locus  + '.map')
    try:
        os.remove(ped)
        os.remove(map_f)
    except OSError:
        logging.warning('Could not remove Plink INPUT files, have they already been removed')
        pass 

def plink_to_ld_matrix(locus ,output_directory , remove_plink_files=False):
    """
        Converts the plink file to a LD matrix for use downstream in paintor.
    """
    output_file = locus + '.LD'
    command = __PLINK_TO_LD_MATRIX__.format(locus,locus)
    # TODO: Fix this workaround, which has to change directory because of a limitation in plink
    try:
        os.chdir(output_directory)
    except OSError:
        logging.error("Could not change directory")
        sys.exit(OS_ERROR)
    run_command(command)
    os.rename(locus + '.ld', locus  + '.matrix')
    _add_dimensions_to_file(locus + '.matrix')
    if remove_plink_files:
        _remove_plink_files(output_directory, locus)
    try:
        os.chdir('../')
    except OSError:
        logging.error("Could not change directory")
        sys.exit(OS_ERROR)
    try:
        os.remove(os.path.join(output_directory,locus +'.log'))
        os.remove(os.path.join(output_directory,locus + '.nosex'))
    except OSError:
        logging.warning('Could not remove Plink INPUT files, have they already been removed')
        pass 


