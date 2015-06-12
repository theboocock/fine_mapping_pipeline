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
#


import os
__ONE_THOUSAND_GENOMES_SAMPLE_MAP__=os.path.join(os.path.dirname(__file__), "../../1000genomes_super_pop.txt")
import logging


def _load_one_thousand_genomes_sample_dict():
    """
        Load the 1000 thousand genomes dataset

        The format is as follows.

        <SAMPLE NAME> <POPULATION>
        EUR NA12839
    """
    one_thousand_genomes_dict = {}
    with open(__ONE_THOUSAND_GENOMES_SAMPLE_MAP__) as samples:
        for sample_l in samples:
            sample_l = sample_l.strip()
            s_pop = sample_l.split('\t')[1]
            sample_name = sample_l.split('\t')[0]
            try:
                one_thousand_genomes_dict[s_pop].append(sample_name)
            except KeyError:
                one_thousand_genomes_dict[s_pop] = [sample_name]
    return one_thousand_genomes_dict

def _get_samples_indices(samples, super_population):
    """
        Obtain the indices to keep from each line of the VCF.

    """
    onekg_dict = _load_one_thousand_genomes_sample_dict()
    super_pop_list = onekg_dict[super_population]
    indices = []
    for i, sample in enumerate(samples):
        if sample in super_pop_list:
            indices.append(i)
    # Let's make sure we return all the indices to keep
    # Need to get columns 1:9
    indices = [ i + 9 for i in indices]
    indices = range(0,9) + indices
    return indices
         

def extract_population_from_1000_genomes(vcf, super_population="EUR", biallelic_only=True, min_maf=0.02):
    """
        Extract a population from a VCF file.

        Function also removes any tri-allelic SNPs
    """
    vcf_temp = ''
    logging.info("Extracting {0} population from VCF".format(super_population))
    for line in vcf.splitlines():
        if "#" in line:
            if "#CHROM" in line:
                samples = line.split('\t')[9:len(line.split('\t'))]
                sample_indices = _get_samples_indices(samples, super_population)
                vcf_temp += '\t'.join([item for i ,item in enumerate(line.split('\t')) if i in sample_indices])  + '\n' 
            else:
                vcf_temp += line + '\n' 
        else:
            vcf_temp_l = None
            if biallelic_only:
                alt = line.split('\t')[4]
                if alt in ['A', 'G', 'C', 'T']:
                    vcf_temp_l = [item for i, item in enumerate(line.split('\t')) if i in sample_indices]
            else:
                vcf_temp_l = [item for i, item in enumerate(line.split('\t')) if i in sample_indices]
            if vcf_temp_l is not None:
                num_aa = len([item for item in vcf_temp_l[9:] if item == '0|0'])
                num_ab = len([item for item in vcf_temp_l[9:] if item == '0|1'])
                num_ab2 = len([item for item in vcf_temp_l[9:] if item == '1|0'])
                num_ab += num_ab2
                num_bb = len([item for item in vcf_temp_l[9:] if item == '1|1'])
                if num_aa == 0 and num_ab == 0:
                    continue
                elif num_ab == 0 and num_bb == 0:
                    continue
                else:
                    numa = num_aa + num_ab 
                    numb = num_ab + num_bb 
                    total_alleles = num_aa + num_ab + num_bb
                    if numa > numb:
                        maf = numa/float(total_alleles)
                    else:
                        maf = numb/float(total_alleles)
                    if maf > min_maf:
                        vcf_temp += '\t'.join(vcf_temp_l) + '\n'

    return vcf_temp

if __name__ == "__main__":
    import doctest
    docetst.testmod()
