import logging

logger = logging.getLogger(__name__)

def chrom_to_number(chrom):
    """
        Converts chromosome based position, such as chr1,chr2
        etc,
        to 1,2
    """
    return chrom.split('chr')[1]
