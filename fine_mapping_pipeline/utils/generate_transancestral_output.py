"""
    Functions to process each locus and transform the LD matrices to fill in 0's where PAINTOR excepts them.

    Also generates overall Z-score file, which is the outer product of the two original files. 

    @author James Boocock
    @date May/25/2016

    Caveat sometimes we will have SNPs with different allele from two different populations at the same position, which were below the MAF threshold in
    one population, but above in another. But, really how many times will this occur. If Z-score files had ref alt this would minimize that. 

    Algorithm:

        1. Create overall SNP list, outer join zscore files
        2. Outer join LD using the same ordering. Filling in the correct zeroes in the columns that need them.
        3. Output new matrices and zscores ready for paintor.
        4. Annotate final Z-score files, which are the merge of the other Z-score files. We will need to ensure we re-run the annotations. 
           

    Should be possible from the original datasets.

"""

import os
import logging
import numpy as np
from orderedset import OrderedSet

class TransAncOutputRow(object):
    """
        Represents transancestral output, needed class so that items 
        could be addede when space was available.
    """

    def __init__(self, chrom, pos, rsid, zscore):
        self._chrom = chrom
        self._pos = str(pos)
        self._rsid = rsid
        self._zscore = zscore 
    
    @property
    def chrom(self):
        return self._chrom
    @property
    def pos(self):
        return self._pos
    @property
    def rsid(self):
        return self._rsid
    @property
    def zscore(self):
        return self._zscore

def generate_transancestral_output(loci, populations, output_directory, test_suffix=""):
    """ 
        Create transancestral matrices for PAINTOR, need to expand the 
        smaller matrices with zeroes so that PAINTOR does not crash ....
   
        @author James Boocock
    """
    np.set_printoptions(suppress=True)
    logging.info("Generating Transancestral Matrices")
    for locus in loci: 
        snp_list_overall = []
        population_lists = {}
        output_rows_population = {}
        for population in populations:
            output_rows_population[population] = {}
            population_lists[population] = []
            with open(os.path.join(output_directory, locus + '.' + population)) as input_snps:
                for line in input_snps:
                    l_s = line.split()
                    chrom = l_s[0]
                    position = int(l_s[1])
                    rsid = l_s[2]
                    z_score = l_s[3]
                    output_rows_population[population][rsid] = TransAncOutputRow(chrom, position, rsid,z_score)
                    population_lists[population].append(rsid)
                    snp_list_overall.append((rsid,position))
        snp_list_overall.sort(key=lambda x: x[1])
        rsids = [x[0] for x in snp_list_overall]
        o_set = OrderedSet(rsids)
        no_snps = len(o_set)
        logging.info("Beginning to create merge ZScore output")
        with open(os.path.join(output_directory, locus  +  test_suffix), "w") as out_zscore:
            out_zscore.write("CHR POS SNP_ID" + ' ' + ' '. join(populations) + "\n")
            for rsid in o_set:
                zscores = []
                for population in populations:
                    try:
                        output_row =  (output_rows_population[population][rsid])
                        chrom = output_row.chrom
                        pos = output_row.pos
                        zscores.append(output_row.zscore)
                    except KeyError:
                        zscores.append("NA")
                out_zscore.write(chrom + " " + pos + " " + rsid + " " + " ".join(zscores) + "\n")
        logging.info("Completed creation of Zscore output")
        logging.info("Started generation of new LD matrices")
        for population in populations:
            output_ld = os.path.join(output_directory, locus + '.LD.' + population)
            plink_ld_matrix= np.loadtxt(output_ld)
            ld_matrix = np.zeros([no_snps, no_snps])
            final_rows = list(o_set)
            rsids = population_lists[population]
            numpy_index_one = 0
            for i in range(no_snps):
                numpy_index_two = 0
                for j in range(no_snps):
                    if o_set[j] in rsids and o_set[i] in rsids:
                        ld_matrix[i,j] = plink_ld_matrix[numpy_index_one,numpy_index_two]
                        numpy_index_two += 1
                    elif i == j:
                        ld_matrix[i, j] = 1 
                if (o_set[i] in rsids):
                    numpy_index_one += 1
            os.rename(output_ld, output_ld + '.old')
            np.savetxt(output_ld + test_suffix, ld_matrix)
        logging.info("Successfully created LD score matrices")

