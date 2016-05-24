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
import logging
import tempfile

from fine_mapping_pipeline.utils.shell import run_command_return_output, run_command
from fine_mapping_pipeline.config import __1000_genomes_sample_map__

def _load_one_thousand_genomes_sample_dict():
    """
        Load the 1000 thousand genomes dataset

        The format is as follows.

        <SAMPLE NAME> <POPULATION>
        EUR NA12839
    """
    one_thousand_genomes_dict = {}
    with open(__1000_genomes_sample_map__) as samples:
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

def _get_sample_list(super_population):
    """
        Loads 1KG dictionary and extract the super population that we are working with. 
    """

    onekg_dict = _load_one_thousand_genomes_sample_dict()
    super_pop_list = onekg_dict[super_population]
    sup_pop_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    for sample in super_pop_list:
        sup_pop_file.write(sample + "\n")
    return sup_pop_file.name

__BCFTOOLS__COMMAND__="""
    bcftools view --force-samples -m2 -M2 -S {0} {1} | bcftools filter -i "MAF > 0.01" 
"""

def _get_cnv_alternate():
    """
        Get an alternate allele when you have a CNV
    """
    return 'A', 'T'


def extract_population_from_1000_genomes(vcf, super_population="EUR", biallelic_only=True, min_maf=0.05):
    """
        Extract a population from a VCF file.

        Function also removes any tri-allelic SNPs
    """
    logging.info("Extracting {0} population from VCF".format(super_population))
    sample_list_file = _get_sample_list(super_population)
    bcftools_command = __BCFTOOLS__COMMAND__.format(sample_list_file, vcf)
    logging.info(bcftools_command)
    #sys.exit(1)
    output_vcf = run_command_return_output(bcftools_command, shell=True)
    last_pos = -1
    vcf_temp = ""
    for i, line in enumerate(output_vcf.splitlines()):
        if i % 1000 ==0 :
            logging.info("Processed {0} lines".format(i))
        if "#" in line:
            vcf_temp += line + '\n' 
        else:
            line_split = line.split("\t")
            rsid = line_split[2]
            position = line_split[1]
            pos2 = int(position)
            if pos2 == last_pos:
                continue
            last_pos = pos2
            reference = line_split[3]
            alternate = line_split[4]
            # Indels and Copy-number
            if "CN" in alternate:
                # Must be a CNV
                reference, alternate = _get_cnv_alternate()
                rsid = "rsCNV" + position
            elif len(reference) > 1 or len(alternate) > 1:
                reference, alternate = _get_cnv_alternate()
                rsid = "rsINDEL" + position
            if rsid == '.':
                rsid = "rs" + chrom + ":" + position + "_" + reference + "/" + alternate
            line_split[3] = reference
            line_split[4] = alternate
            line_split[2] = rsid
            vcf_temp += "\t".join(line_split) + "\n"
    #for i, line in enumerate(vcf.splitlines()):
    #    if (i %1000 == 0):
    #        logging.info("Processed {0} lines from the VCF file".format(i))
    #    if "#" in line:

    #        if "#CHROM" in line:
    #            samples = line.split('\t')[9:len(line.split('\t'))]
    #            sample_indices = _get_samples_indices(samples, super_population)
    #            vcf_temp += '\t'.join([item for i ,item in enumerate(line.split('\t')) if i in sample_indices])  + '\n' 
    #        else:
    #            vcf_temp += line + '\n' 
    #    else:
    #        vcf_temp_l = None
    #        if biallelic_only:
    #            alt = line.split('\t')[4]
    #            if "," not in alt:
    #                vcf_temp_l = [item for i, item in enumerate(line.split('\t')) if i in sample_indices]
    #        else:
    #            vcf_temp_l = [item for i, item in enumerate(line.split('\t')) if i in sample_indices]
    #        if vcf_temp_l is not None:
    #            num_aa = len([item for item in vcf_temp_l[9:] if item == '0|0'])
    #            num_ab = len([item for item in vcf_temp_l[9:] if item == '0|1'])
    #            num_ab2 = len([item for item in vcf_temp_l[9:] if item == '1|0'])
    #            num_ab += num_ab2
    #            num_bb = len([item for item in vcf_temp_l[9:] if item == '1|1'])
    #            if num_aa == 0 and num_ab == 0:
    #                continue
    #            elif num_ab == 0 and num_bb == 0:
    #                continue
    #            else:
    #                numa = num_aa + num_ab 
    #                numb = num_ab + num_bb 
    #                total_alleles = num_aa + num_ab + num_bb
    #                if numa > numb:
    #                    maf = numa/float(total_alleles)
    #                else:
    #                    maf = numb/float(total_alleles)
    #                if maf > min_maf:
    #                    vcf_temp += '\t'.join(vcf_temp_l) + '\n'

    return vcf_temp

if __name__ == "__main__":
    import doctest
    docetst.testmod()
