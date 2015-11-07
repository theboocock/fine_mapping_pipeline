# Repair correlatino matrices.

Basically if a SNP is a proxy to another SNP it should be removed.

This should be easy enough, scan along the rows and look for 1's or negative 1's
in anything other than the diagonal.

Ideally this should be incorporated into the main dataset

# Zscore to pvalue.

Takes a list of Zscores from a region with RSIDS followed by Zscore and converts them to Pvalues.


