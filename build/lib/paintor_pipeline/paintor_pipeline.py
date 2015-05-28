#!/usr/bin/env python
#
import argparse
import os
from pyrallel import *
from paintor_pipeline.ucsc_utilities import * 

onekg_index = os.path.join(os.path.abspath(__file__), "file_index/vcf_files_1kg.txt")

class snpPos:
    def __init__(self,chrom, pos):
        self.chrom=chrom
        self.pos=pos
    @property
    def chrom(self):
        return self.chrom
    @property
    def pos(self):
        return self.pos

def get_snp_list(snp_list):
    # Check if the RSID is working.
    snps_position=[]
    with open(snp_list) as f:
        for line in f:
            line = line.strip()
            if "rs" in line: 
                (chrom, pos) = get_position_ucsc(line,build)
            else:
                chrom = line.split('\t')[0]
                pos = line.split('\t')[1]
            snps_position.append(snpPos(chrom,pos))      



def main():
    """
        Process SNP for a Paintor run. 
    """
    parser = argparse.ArgumentParser(description="Runs paintor across a list of SNPS")
    parser.add_argument('-s','--snp-list', dest='snp_list_file', help='Contains a list of chrom:position files or rsids')
    parser.add_argument('-z','--z-score-dir', dest='z-score-dir', help='Z-score directory')
    parser.add_argument('-b','--build',dest='genome_build', help='Genome build (HG19)')
    args = parser.parse_args()
    get_snp_list(args.snp_list_file)

if __name__ == "__main__":
    main()

