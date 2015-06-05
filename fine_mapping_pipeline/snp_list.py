from paintor_pipeline.ucsc import snp_utilities 
from collections import Iterator

class Snp:
    """
        Class reperesents a SNP with a chromosome and position for use in downstream
        analyses.
    """

    def __init__(self, chrom, pos,rs_id):
        self.chrom = chrom
        self.pos = int(pos)
        self.rsid = rs_id

    @property 
    def chrom(self):
        return self.chrom
    @property 
    def pos(self):
        return self.pos
    @property
    def rsid(self):
        return self.rsid

    def __str__(self):
        return self.chrom+":"+str(self.pos)

class SnpList(Iterator):

    def __init__(self,snp_list, build):
        """
            Create SNP list from file
        """
        self.build = build
        self.snp_list = self._create_snp_list(snp_list)

    def _create_snp_list(self, snp_list):
        """ 
            Process SNP list file line by line
        """
        list_of_snp_objs = []
        with open(snp_list) as snps:
            for line in snps:
                rs_maybe = line.split()[0]
                if "rs" in rs_maybe:
                    (chrom, pos) =  snp_utilities.get_chrom_pos_from_rsid(rs_maybe, self.build)
                else:
                    chrom = line.split()[0]
                    pos = line.split()[1]
                list_of_snp_objs.append(Snp(chrom,pos,rs_maybe))
        return list_of_snp_objs

    def __str__(self):
        """
            To string method for a snp list
        """
        output_l =[] 
        for snp in self.snp_list:
             output_l.append(str(snp))
        return '\n'.join(output_l)

    def next(self):
        """
            Next method for a SNP list.
        """
        if not self.snp_list:
            raise StopIteration
        return self.snp_list.pop()
