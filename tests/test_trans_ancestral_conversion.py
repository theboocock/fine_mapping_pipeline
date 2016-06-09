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

import os

from fine_mapping_pipeline.utils.generate_transancestral_output import generate_transancestral_output

def test_generate_transancestral_runs():
    """
        Test that we can generate transancestral output
    """
    transancestral_dir = "tests/trans_ances"
    populations = ["EUR","EAS"]
    loci = ['1_4:88800000-88900000']
    test_suffix = ".tmp"
    generate_transancestral_output(loci, populations, output_directory=transancestral_dir, test_suffix=test_suffix) 
    with open(os.path.join(transancestral_dir, loci[0] + '.LD.' + populations[0] +".tmp")) as f:
        for i, line in enumerate(f):
            pass
        assert i == 520


def test_generate_transancestral_interpolate_matrix():
    """
        Test generate transancestral interpolate matrix function using a simple example.
    """
    transancestral_dir = "tests/trans_ances2"
    populations = ["EUR","EAS"]
    loci = ['1_4:88800000-88900000']
    generate_transancestral_output(loci, populations, output_directory=transancestral_dir, test_suffix=test_suffix) 
    with open(os.path.join(transancestral_dir, loci[0] + '.LD.' + populations[0] +".tmp")) as f:
        for i, line in enumerate(f):
            if (i ==0):
                a_line = line.split()
                a_line = [float(a) for a in a_line]
                assert a_line[0] == 1.0
                assert a_line[1] == 1.0
                assert a_line[2] == 0.0
                assert a_line[3] == 0.5
                assert a_line[4] == 0.0
        assert i == 4
    with open(os.path.join(transancestral_dir, loci[0] + '.LD.' + populations[1] +".tmp")) as f:
        for i, line in enumerate(f):
            if (i ==0):
                a_line = line.split()
                a_line = [float(a) for a in a_line]
                assert a_line[0] == 1.0
                assert a_line[1] == 1.0
                assert a_line[2] == 0.0
                assert a_line[3] == 0.0
                assert a_line[4] == 0.5
        assert i == 4
