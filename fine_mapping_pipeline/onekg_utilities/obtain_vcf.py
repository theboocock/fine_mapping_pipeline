__1000_genomes_template__="""ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"""

import subprocess
import shlex
import os
import sys

import logging 
logger = logging.getLogger(__name__)
from fine_mapping_pipeline.expections import *
from fine_mapping_pipeline.utils.shell import *

def get_vcf_file(snp, flanking_region):
    """
        Obtain a SNP and a flanking region
    """
    
    vcf = __1000_genomes_template__.format(snp.chrom)
    chrom = snp.chrom
    pos = snp.pos 
    start = snp.pos - flanking_region
    end = snp.pos + flanking_region
    command = 'tabix -h ' + vcf + ' ' + chrom + ':' + str(start) + '-' + str(end)
    vcf_file_in_memory = run_command_return_output(command)
    return vcf_file_in_memory
