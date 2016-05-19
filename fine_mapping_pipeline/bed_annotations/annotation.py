# Copyright (c) 2016 Boocock James <james.boocock@otago.ac.nz>
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
import numpy as np
from bedtools import IntervalFile
 
class AnnotateLociMatrix(object):
    """
        Run the entirety of the annotations.

        @date May 17 2016
    """

    def __init__(self, no_annotiations, no_snps):
        self._zeroes_=  np.zeroes([no_snps,no_annotiationsnnotations])
        self._header = []
        

def _bed_from_vcf(vcf):
    """
        Bed file from intersection
    """
    with open(vcf) as vcf_input_file:

def _get_line_number(vcf):
    """ 
        Get number of lines in vcf file 
    """
    with open(vcf) as vcf_input_file:
        i =  0 
        for line in vcf_input_file:
            if "#" in line:
                i += 1
    return i

def generate_bed_file_annotations(bed_directory, output_directory, vcf):
    """
        Generates the annotation file for every bed file in the bed_directory folder
    """
    
    # Loop over the bed files in the bed directory.
    tmp_bed = _bed_from_vcf(vcf)
    no_snps = _get_line_number(vcf)
    bed_file_list = glob.glob(os.path.join(args.bed_directory, "*.bed"))
    a_matrix= AnnotateMatrix(no_snps)
    for beds in bed_file_list: 

