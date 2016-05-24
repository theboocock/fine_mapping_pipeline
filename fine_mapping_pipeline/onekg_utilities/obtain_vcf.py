
import subprocess
import shlex
import os
import sys

import logging

import tempfile

logger = logging.getLogger(__name__)
from fine_mapping_pipeline.expections import *
from fine_mapping_pipeline.utils.shell import *
from fine_mapping_pipeline.config import __1000_genomes_template__


def remove_tbi_files():
    """ 
        Should probably use this function to remove any temp files that create downstream problems.
    """
def get_vcf_file(snp, flanking_region=None, string_region=None):
    """
        Obtain a SNP and a flanking region
    """
    vcf = __1000_genomes_template__.format(snp.chrom)
    chrom = snp.chrom
    pos = snp.pos
    if flanking_region is not None:
        start = pos - flanking_region
        end = pos + flanking_region
        command = 'tabix -h ' + vcf + ' ' + chrom + ':' + str(start) + '-' + str(end)
    elif string_region is not None:
        rsid = snp.rsid
        chrom = string_region.split(":")[0]
        positions = string_region.split(":")[1]
        start = positions.split("-")[0]
        end = positions.split("-")[1]
        command = 'tabix -h ' + vcf + ' ' + chrom + ':' + str(start) + '-' + str(end)
    tf = tempfile.NamedTemporaryFile(delete=False, mode="w")
    run_command(command, stdout=tf)
    logging.info(tf.name)
    return tf.name 
