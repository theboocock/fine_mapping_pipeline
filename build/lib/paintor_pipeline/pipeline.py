import argparse
import paintor_pipeline.ucsc

from paintor_pipeline.expections import error_codes

from paintor_pipeline.snp_list import SnpList
from paintor_pipeline.onekg_utilities.obtain_vcf import get_vcf_file 

def paintor_run(args):
    """
        Parses arguments from the paintor sub-command and processes data for use in Paintor.

    """
    z_score_dir = args.z_score_dir
    try:
        flanking_region = int(args.flanking_region)
    except ValueError:
        logger.error('Flanking region argument needs to be an integer')
        sys.exit(COMMAND_LINE_ERROR)
    build = args.build
    
    # Create the SNPList
    snp_list = SnpList(args.snp_list,build)
    for snp in snp_list:
        vcf = get_vcf_file(snp, flanking_region)


def main():
    """
        Creates and runs a fine mapping analysis.

    """
    parser = argparse.ArgumentParser(description="Processes SNP based data and performs various fine mapping tasks")
    subparsers = parser.add_subparsers(help='Sub-command help')
    # Paintor parser 
    paintor_parser = subparsers.add_parser('paintor',help='Paintor Help')
    paintor_parser.add_argument('-s','--snp-list', dest='snp_list', help='SNP List file rsids or bed formatted')
    paintor_parser.add_argument('-z','--z-score-dir', dest='z_score_dir', help='File containing Z-scores for an entire chromosome of SNPs')
    paintor_parser.add_argument('-f','--flanking-region', dest='flanking_region', help='Flanking region')
    paintor_parser.add_argument('-b','--build', dest='build', help='Genome build',default='hg19')
    paintor_parser.set_defaults(func=paintor_run)
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
