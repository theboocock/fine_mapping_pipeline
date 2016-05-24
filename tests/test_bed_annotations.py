# Copyright (c) 2016 Boocock James <james.boocock@otago.ac.nz>
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

# Test Bed annotations.

from fine_mapping_pipeline.bed_annotations.annotation import generate_bed_file_annotations

def test_bed_annotiations(tmpdir):
    """
        Tests bed annotation of the Zscore files from IMPG 
    """
    bed_directory= "tests/bed_files/"
    vcf_input_files = ["tests/vcfs/small.vcf", "tests/vcfs/test.vcf"]
    loci = ["small", "test"]
    population = "EUR"
    output_directory = str(tmpdir)
    generate_bed_file_annotations(bed_directory=bed_directory, output_directory=output_directory, vcfs=vcf_input_files, loci=loci, population=population)
    assert 1== 2
