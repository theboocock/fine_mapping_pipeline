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


import argparse
import paintor_pipeline.ucsc
## Import datetime to make dated directory
import time, datetime, os

from os import listdir

from paintor_pipeline.expections import error_codes

from paintor_pipeline.snp_list import SnpList
from paintor_pipeline.onekg_utilities.obtain_vcf import get_vcf_file 
from paintor_pipeline.gemini.create import create_gemini_database  
from paintor_pipeline.gemini.annotation import encode_annotation 

def get_relevant_zscore(chrom, directory):
    """
        Scans the directory specified on the command-line by the user and checks for the specified chromosome.

    >>> get_relevant_zscore('6', t
    """
    zscore_files = [ f for f in listdir(directory) ]
    chrom = 'chr' + chrom
    z_score_file = [ f for f in zscore_files if chrom in f ][0]
    return os.path.join(directory, z_score_file)

def create_pos_hash_table(zscore_file):
    """
        Create position hash table

        @param pos_file the Z scores file containing the positions from impg

        @return pos_hash a hash table that contains the line from impG

    """
    pos_hash = {}
    with open(zscore_file) as p:
        for i, line in enumerate(p):
            line = line.strip()
            if i != 0:
                pos_hash[int(line.split('\t')[1])] = line.split('\t')[3]
    return (pos_hash)

def generate_zscore_and_vcf_output(output_directory, 
                             zscore_hash, 
                             vcf,
                             locus):
    """
        Extract vcf regions that have overlap with Impg data.

    """
    # output vcf is a temporary file that will be used downstream of this dataset.
    output_vcf = os.path.join(output_directory, locus+'.vcf') 
    output_zscore = os.path.join(output_directory, locus)
    with open(output_vcf, 'w') as out_vcf:
        with open(output_zscore, 'w') as out_zscore:
            for line in vcf.splitlines():
                if "#" in line:
                    out_vcf.write(line)
                else:
                    pos = int(line.split('\t')[1])
                    t_chrom = line.split('\t')[0]
                    if pos in zscore_hash.keys():
                        out_vcf.write(line)
                        out_zscore.write(zscore_hash[pos]+'\n')

def prepare_run(args):
    """
        Parses arguments from the paintor sub-command and processes data for use in Paintor.

    """
    if args.output_directory is None:
        todaystr = 'paintor_run' + datetime.date.today().isoformat()
        output_directory = todaystr
        try:
            os.mkdir(todaystr)
        except OSError:
            pass
    else:
        output_directory = args.output_directory
    z_score_dir = args.z_score_dir
    try:
        flanking_region = int(args.flanking_region)
    except ValueError:
        logger.error('Flanking region argument needs to be an integer')
        sys.exit(COMMAND_LINE_ERROR)
    build = args.build
    # Create the SNPList
    snp_list = SnpList(args.snp_list, build)
    # Locus to process 
    loci = []
    for snp in snp_list:
        locus = snp.rsid
        loci.append(locus)
        vcf = get_vcf_file(snp, flanking_region)
        z_score_file =  get_relevant_zscore(snp.chrom, z_score_dir)
        pos_list_zscore = create_pos_hash_table(z_score_file) 
        generate_zscore_and_vcf_output(output_directory=output_directory, zscore_hash=pos_list_zscore, vcf=vcf, locus=locus) 
        gemini_database = create_gemini_database(vcf)
        create_gemini_annotations(gemini_database, locus, output_directory)  

def main():
    """
        Creates and runs a fine mapping analysis.

    """
    parser = argparse.ArgumentParser(description="Processes SNP based data and performs various fine mapping tasks")
    subparsers = parser.add_subparsers(help='Sub-command help')
    # Prepare parser 
    prepare_parser = subparsers.add_parser('prepare',help='Prepare for paintor run')
    prepare_parser.add_argument('-s','--snp-list', dest='snp_list', help='SNP List file rsids or bed formatted')
    prepare_parser.add_argument('-z','--z-score-dir', dest='z_score_dir', help='File containing Z-scores for an entire chromosome of SNPs')
    prepare_parser.add_argument('-f','--flanking-region', dest='flanking_region', help='Flanking region')
    prepare_parser.add_argument('-b','--build', dest='build', help='Genome build',default='hg19')
    prepare_parser.add_argument('-o','--output', dest='output_directory', help='output directory empty or non-existent directory for dumping files to be used in a paintor run')
    prepare_parser.add_argument('-p','--population', dest='population', help='1kg population to calculate LD from', default='EUR')
    prepare_parser.add_argument('-m','--maf',dest='maf', help='MAF filtering for 1000 genomes VCF file', default=0.01)
    prepare_parser.set_defaults(func=prepare_run)
    
    #Paintor run parser
    paintor_parser = subparsers.add_parser('paintor', help='Sub command runs paintor following file preparation')
    paintor_parser.add_argument('-d','--input_directory',dest='input_directory', help="Directory files were prepared in after running the" 
                                "prepare command")
    paintor_parser.add_argument('-a','-auto-annotations', help='Run paintor multiple times to select the annotations to use')

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
