
import subprocess
import shlex
import os
import sys

import logging 
logger = logging.getLogger(__name__)
from fine_mapping_pipeline.expections import *
from fine_mapping_pipeline.utils.shell import *
from fine_mapping_pipeline.config import __1000_genomes_template__

def remove_tbi_files():
    """ 
        Should probably use this function to remove any temp files that create downstream problems.
    """
def get_vcf_file(snp, flanking_region):
    """
        Obtain a SNP and a flanking region
    """
    vcf = __1000_genomes_template__.format(snp.chrom)
    chrom = snp.chrom
    pos = snp.pos
    start = pos - flanking_region
    end = pos + flanking_region
    command = 'tabix -h ' + vcf + ' ' + chrom + ':' + str(start) + '-' + str(end)
    vcf_file_in_memory = run_command_return_output(command)
    return vcf_file_in_memory
