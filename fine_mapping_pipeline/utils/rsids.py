"""
    Misc RSID processing to keep it consistent.

    Basically parse RSIDs and place them.
"""

def normalise_rsids(rsids, pos, reference, alternate):
    """
        Normalise an rsid
    """
    if "CN" in alternate:
        # Must be a CNV
        reference, alternate = _get_cnv_alternate()
        rsid = "rsCNV" + position
    elif len(reference) > 1 or len(alternate) > 1:
        reference, alternate = _get_cnv_alternate()
        rsid = "rsINDEL" + position
    if rsid == '.':
        rsid = "rs" + chrom + ":" + position + "_" + reference + "/" + alternate
    rsid = line_split[2]
    position = line_split[1]
