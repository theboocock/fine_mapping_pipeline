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

import fnmatch
import fine_mapping_pipeline.ucsc
## Import datetime to make dated directory
import time, datetime, os, sys

import pandas as pd 

from fine_mapping_pipeline.expections.error_codes import *

from fine_mapping_pipeline.plink.ld import vcf_to_plink, plink_to_ld_matrix
from fine_mapping_pipeline.snp_list import SnpList, Snp
from fine_mapping_pipeline.onekg_utilities.obtain_vcf import get_vcf_file
from fine_mapping_pipeline.onekg_utilities.vcf_filter import extract_population_from_1000_genomes
from fine_mapping_pipeline.gemini.create import create_gemini_database
from fine_mapping_pipeline.gemini.annotation import generate_and_write_encode_annotations
from fine_mapping_pipeline.bed_annotations.annotation import generate_bed_file_annotations 

from fine_mapping_pipeline.utils.zscores import get_relevant_zscore, create_pos_hash_table, generate_zscore_and_vcf_output
from fine_mapping_pipeline.finemap.paintor import run_paintor
from fine_mapping_pipeline.finemap.caviarbf import run_caviarbf
from fine_mapping_pipeline.utils.generate_transancestral_output import generate_transancestral_output 


def _prepare_output_dir(output_directory):
    if output_directory is None:
        todaystr = 'fine_mapping_run' + datetime.date.today().isoformat()
        output_directory = todaystr
    else:
        output_directory = output_directory
    try:
        os.mkdir(output_directory)
    except OSError:
        pass
    if not output_directory.endswith('/'):
        output_directory += '/'
    return output_directory

def prepare_runs(args):
    """
        Parses arguments from the paintor sub-command and processes data for use in Paintor.

    """
    output_directory = _prepare_output_dir(args.output_directory)
    z_score_dir = args.z_score_dir
    region_list = args.region_list 
    if args.region_list is None:
        try:
            flanking_region = int(args.flanking_region)
        except ValueError:
            logging.error('Flanking region argument needs to be an integer')
            sys.exit(COMMAND_LINE_ERROR)
    build = args.build
    bed_directory = args.bed_directory
    # Create the SNPList
    try:
        min_maf = float(args.maf)
    except:
        logging.error("Min Maf -m or --min-maf needs to be an floating point number")
        sys.exit(COMMAND_LINE_ERROR)
    if args.region_list is not None:
        region_list = {}
        snp_list = []
        with open(args.region_list) as input_file:
            # When using no flaking region SNP must be valid, but it doesn't actually matter what it is, need to ensure that is actually the case.
            for i, line in enumerate(input_file):
                rsid = str(i)+ "_"  + ''.join(line.strip().split("\t"))
                chromosome = line.strip().split(":")[0] 
                snp = Snp(chromosome,"1",rsid)
                snp_list.append(snp)
                region_list[snp.rsid] = line.strip()
    else:
        snp_list = SnpList(args.snp_list, build)
        logging.info(snp_list)
    # Locus to process
    # population_to_extract_vcf
    if not args.annotation_only:
        no_flanking = args.flanking_units
        if no_flanking:
            raise NotImplementedError("Using a number of flanking SNPs instead of a region is not supported")
        populations= args.populations.split(',')
        logging.info("Populations to process: {0}".format(populations))
        loci = []
        gemini_databases = []
        output_vcfs = []
        for snp in snp_list:
            logging.info('Preparing output files for SNP {0}'.format(snp.rsid))
            locus = snp.rsid
            loci.append(locus)
            logging.info("Obtaining VCF file from the 1000 genomes project")
            if region_list is not None:
                vcf = get_vcf_file(snp, string_region=region_list[locus])
            else:    
                vcf = get_vcf_file(snp, flanking_region=flanking_region)
            for population in populations:
                tmp_vcf = extract_population_from_1000_genomes(vcf=vcf, super_population=population)
                z_score_file = get_relevant_zscore(snp.chrom, population, z_score_dir)
                pos_list_zscore = create_pos_hash_table(z_score_file)
                output_vcf = generate_zscore_and_vcf_output(output_directory=output_directory, zscore_hash=pos_list_zscore, vcf=tmp_vcf, locus=locus,population=population, multiply_rsquare=args.multiply_rsquare)
                if bed_directory is None:
                    logging.info("Creating gemini database")
                    # TODO: Fix broxen gemini referenec
                    gemini_databases.append(create_gemini_database(vcf=output_vcf))
                vcf_to_plink(locus, output_directory=output_directory, vcf=output_vcf, population=population)
                plink_to_ld_matrix(locus, output_directory=output_directory, population=population)
        logging.info("Generate transancestrals matrices")
        generate_transancestral_output(loci, populations, output_directory)
        if bed_directory is None:
            logging.info("Generating annotation matrices to be used with Paintor")
            logging.info(gemini_databases)
            generate_and_write_encode_annotations(databases=gemini_databases, output_directory=output_directory, loci=snp_list)
        else:
            logging.info("Annotation using bed files")
            generate_bed_file_annotations(loci=loci, bed_directory=bed_directory, output_directory=output_directory) 
        # So finally we need to fix the LD matrices for inputting into PAINTOR. 

        with open(os.path.join(output_directory, 'input.files'), 'w') as out_f:
            for snp in snp_list:
                out_f.write(snp.rsid +'\n')
        # Remove .tbi files
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, '*.tbi'):
            try:
                os.remove(file)
            except OSError:
                logging.warning("Could not remove a .tbi file from the 1000 genomes tabix run")
    else: 
        loci = []
        for snp in snp_list:
            loci.append(snp.rsid)
        if bed_directory is not None:
            logging.info("Annotation using bed files")
            generate_bed_file_annotations(loci=loci, bed_directory=bed_directory, output_directory=output_directory) 
    logging.info("Finemapping file preparation complete")

