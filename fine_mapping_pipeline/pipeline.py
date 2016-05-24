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

import sys
import logging
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s %(levelname)s %(message)s')
import argparse
from fine_mapping_pipeline.prepare_input.prepare_runs import prepare_runs
from fine_mapping_pipeline.finemap.caviarbf import run_caviarbf_wrap 
from fine_mapping_pipeline.finemap.paintor import run_paintor_wrap
from fine_mapping_pipeline.finemap.finemap import run_finemap_wrap

def main():
    """
        Creates and runs a fine mapping analysis.

    """
    logging.info('Starting the fine mapping pipeline')
    parser = argparse.ArgumentParser(description="Processes SNP based data and performs various fine mapping tasks")
    subparsers = parser.add_subparsers(help='Sub-command help')
    # Impg parser

    # Prepare parser 
    prepare_parser = subparsers.add_parser('prepare',help='Prepare for paintor run')
    prepare_parser.add_argument('-s', '--snp-list', dest='snp_list', help='SNP List file rsids or bed formatted')
    prepare_parser.add_argument('-z', '--z-score-dir', dest='z_score_dir', help='File containing imputed Z-scores for an entire chromosome of SNPs')
    prepare_parser.add_argument('-f', '--flanking-region', dest='flanking_region', help='Flanking region')
    prepare_parser.add_argument('-n', '--number_of_snps', dest='flanking_units',action='store_true',help='Use a number of SNPs either side instead of a region', default=False)
    prepare_parser.add_argument('-g', '--genome-build', dest='build', help='Genome build',default='hg19')
    prepare_parser.add_argument('-o', '--output', dest='output_directory', help='output directory empty or non-existent directory for dumping files to be used in a paintor run')
    prepare_parser.add_argument('-p', '--populations', dest='populations',help='1kg population to calculate LD from', default='EUR')
    prepare_parser.add_argument('-m', '--min-maf', dest='maf', help='Min MAF filtering for 1000 genomes VCF file', default=0.02)
    prepare_parser.add_argument("-b", "--bed-directory", dest="bed_directory", help="Bed file directory (does not utilise gemini)")
    prepare_parser.add_argument('-r', '--region_list', dest='region_list', help='Region list file', required=False) 
    prepare_parser.set_defaults(func=prepare_runs)
    
    #Finemap parser
    paintor_parser = subparsers.add_parser('paintor', help='Run and process paintor\
                                           data following file preparation')
    paintor_parser.add_argument('-i','--input_directory', dest='input_directory', 
                                help="Directory files were prepared in after running the"
                                "prepare command", required=True)
    paintor_parser.add_argument('-a', '-auto-annotations', dest='auto_select_annotations',action='store_true', 
                                help='If using paintor select the annotations.', required=True)
    paintor_parser.add_argument('-o','--output_directory', dest='output_directory', help="Results output dir")
    paintor_parser.add_argument('-c','--causal_snp_number', dest='causal_snp_number', help="Potiential number of casual SNPs",
                                default=3)
    paintor_parser.set_defaults(func=run_paintor_wrap)
    
    caviarbf_parser = subparsers.add_parser('caviarbf', help='Run and process caviarbf\
                                            output following file preparation')
    caviarbf_parser.add_argument('-i','--input_directory', dest='input_directory', 
                                help="Directory files were prepared in after running the\
                                prepare command", required=True)
    caviarbf_parser.add_argument('-o', '--output_directory', dest='output_directory', help="Results output dir")
    caviarbf_parser.add_argument('-c', '--causal_snp_number', dest='causal_snp_number', help='Potential number of causal SNPs',
                                 default=3)
    caviarbf_parser.add_argument('-p', '--prior', dest='prior', help='Caviarbf prior information',
                                 default=0.1281429)
    caviarbf_parser.add_argument('-s', '--sample-size', dest='sample_size', 
                                 help='Sample size Zscores were calculated from', required=True)
    caviarbf_parser.set_defaults(func=run_caviarbf_wrap)

    finemap_parser = subparsers.add_parser("finemap", help='Run and process finemap\
                                            output following file preparation')
    finemap_parser.add_argument('-i','--input_directory', dest='input_directory',
                                help="Directory files were prepare in after running the\
                                prepare command", required=True)
    finemap_parser.add_argument("-o", "--output_directory", dest="output_directory", help="Results output dir")
    finemap_parser.add_argument("-c", "--causal_snp_number", dest="causal_snp_number", help="Potential number of causal SNPs",
                              default=3)
    finemap_parser.add_argument("-n","--n-ind", dest="number_of_individuals", help="Number of individuals")
    finemap_parser.set_defaults(func=run_finemap_wrap)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
