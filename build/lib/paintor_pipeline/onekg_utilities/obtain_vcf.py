__1000_genomes_template__="""
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
"""

import tabix
import subprocess
import shlex
import os

from paintor_pipeline.expections import *

def get_vcf_file(snp, flanking_region):
    """
        Obtain a SNP and a flanking region
    """
    
    vcf = __1000_genomes_template__.format(snp)
    chrom = snp.chrom
    pos = snp.pos 
    start = snp.pos - flanking_region
    end = snp.pos + flanking_region
    command = 'tabix -h ' + vcf + ' ' + chrom + ':' + str(start) + '-' + str(end)
    try:
       command = shlex.split(command)
       vcf_file_in_memory = subprocess.check_output(command)
       os.remove(os.path.basename(vcf) + '.tbi')
    except subprocess.CalledProcessError:
        logger.error("Problem obtaining VCF file from 1000 genomes ftp - check you connection\n") 
        sys.exit(ONEKG_DOWNLOAD_FAILED)
    return vcf_file_in_memory
