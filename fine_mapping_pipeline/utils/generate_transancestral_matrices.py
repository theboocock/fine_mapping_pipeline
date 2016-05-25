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




def generate_transancestral_matrix(locus, populations):
    """ 
        Function process one locus at a time.
        @date May 25 2016
        @author James Boocock
    """
