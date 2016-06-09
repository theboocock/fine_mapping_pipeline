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
import glob
import os
import tempfile
from pybedtools import BedTool 
 
class AnnotateLociMatrix(object):
    """
        Run the entirety of the annotations.

        @date May 17 2016
    """

    def __init__(self, no_annot , no_snps):
        self._zeroes = np.zeros([no_snps, no_annot],dtype=np.int32)
        self._header = []
        self._annot_idx = 0

    def add_annotation(self, column, annot_name):
        """
            Add annotations to the final analysis file.
        """
        annot_name = os.path.basename(annot_name)
        self._header.append(annot_name)
        self._zeroes[:,self._annot_idx] = column
        self._annot_idx += 1

    def write_annotations(self, output_file):
        """
            Write the annotation matrices
        """
        logging.info(self._header)
        np.savetxt(output_file, self._zeroes, header=" ".join(self._header),fmt='%i',comments='')

#def _bed_from_vcf(vcf):
#    """
#        Bed file from intersection
#    """
#    bed_lines = []
#    rsids = []
#    with open(vcf) as vcf_input_file:
#        for line in vcf_input_file:
#            if "#" not in line:
#                l_s = line.split("\t")   
#                chrom = l_s[0] 
#                pos = l_s[1]
#                rsid = l_s[2]
#                rsids.append(rsid)
#                bed_line = "chr"+l_s[0] + "\t" + str(int(pos) -1 )+ "\t" + str(pos)+ "\t" + rsid +"\n"
#                bed_lines.append(bed_line)
#    return bed_lines, rsids


def _bed_from_zscore(zscore):
    """
        Convert zscore file to bedfile for intersection processing. 
    """
    bed_lines = [] 
    rsids = [] 
    with open(zscore) as zscore_iter:
        for i, line in enumerate(zscore_iter):
            if i > 0:
                l_s = line.split()
                chrom = l_s[0]
                pos = l_s[1]
                rsid = l_s[2]
                rsids.append(rsid)
                bed_line = "chr" + chrom + "\t" + str(int(pos) - 1) + "\t" + str(pos) + "\t" + rsid + "\n"
                bed_lines.append(bed_line)
    return bed_lines, rsids


def _get_line_number(vcf):
    """ 
        Get number of lines in vcf file 
    """
    with open(vcf) as vcf_input_file:
        i =  -1
        for line in vcf_input_file:
                i += 1
    return i

def generate_bed_file_annotations(bed_directory, output_directory, loci):
    """
        Generates the annotation file for every bed file in the bed_directory folder
    """
    
    # Loop over the bed files in the bed directory.
    bed_file_list = glob.glob(os.path.join(bed_directory, "*.bed"))
    logging.info("Start to generate BED file annotations")
    logging.info("Writing annotation to: {0}/".format(output_directory))
    for locus in loci:
        zscore = os.path.join(output_directory, locus) 
        bed_lines, rsids = _bed_from_zscore(zscore)
        tmp_bed = open("tmp.bed","w").writelines(bed_lines)
        snps = BedTool("tmp.bed")
        no_snps = _get_line_number(zscore)
        a_matrix= AnnotateLociMatrix(len(bed_file_list), no_snps)
        logging.info("Annotating locus: {0}, using VCF file {1}".format(locus, zscore))
        for beds in bed_file_list:
            test_annotation = BedTool(beds)
            inter = snps.intersect(test_annotation)
            idxs = []
            for inte in inter:
                idxs.append(rsids.index(inte.name))
            zeroes = np.zeros(len(rsids))
            for idx in idxs:
                zeroes[idx] = 1
            a_matrix.add_annotation(zeroes, beds)
        annotations_file = os.path.join(output_directory, locus + ".annotations")
        logging.info("Writing annotation matrix to: {0}".format(annotations_file))
        a_matrix.write_annotations(annotations_file)
        os.remove("tmp.bed")
