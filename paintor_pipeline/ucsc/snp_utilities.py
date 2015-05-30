# Query ucsc for information.

import subprocess
import logging
import shlex 
from paintor_pipeline.expections.error_codes import *
from paintor_pipeline.ucsc.utils import *
from paintor_pipeline.utils.shell import *


logger = logging.getLogger(__name__)


## Here chromEnd equals the position for that SNP.
__GET_RSID_QUERY__="""
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -NA {0} \
    -e "select chrom, chromEnd from snp141 where name='{1}'"  
"""

def get_chrom_pos_from_rsid(rsid, build):
    """

    Extracts a SNP and position from UCSC snp141 database

    >>> get_chrom_pos_from_rsid("rs12498742","hg19")
    ('4', '9944052')
    """
    query = __GET_RSID_QUERY__.format(build, rsid)
    output = run_command_return_output(query, error=UCSC_QUERY_FAILED)
    output = output.strip() 
    chrom = output.split('\t')[0]
    chrom = chrom_to_number(chrom) 
    position = output.split('\t')[1]
    return chrom, position    



if __name__ == "__main__":
    import doctest
    doctest.testmod()
