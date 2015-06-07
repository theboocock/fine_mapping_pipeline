#!/usr/bin/env bash

fine_mapping_pipeline prepare -s tests/test_data/snp_list/maf_snps.txt -z tests/test_data/zscores/ -f 20000 -b hg19 -o big_test

fine_mapping_pipeline  paintor -i big_test -o big_test/paintor -a
fine_mapping_pipeline  caviarbf -i big_test -o big_test/caviarbf 
