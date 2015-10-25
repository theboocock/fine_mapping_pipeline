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

from fine_mapping_pipeline.finemap.finemap import run_finemap, remove_surrogates, _write_matrix, _write_zscores 
import logging 
logging.basicConfig(level=logging.INFO)

def test_remove_surrogate(tmpdir):
    input_matrix = 'tests/finemap_data/test.matrix'
    input_zscore = 'tests/finemap_data/test.Z'
    surrogates_out = 'tests/finemap_data/out.surro'
    (matrix, zscores) = remove_surrogates(input_matrix,input_zscore, surrogates_out)
    _write_matrix(matrix, "tests/finemap_data/out.matrix")
    _write_zscores(zscores, "tests/finemap_data/out.zscores")
    assert 1 == 2

