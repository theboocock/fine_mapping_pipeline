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

import logging
import os

__PAINTOR_TEMPLATE__="""
    PAINTOR -o {0} -d {1} -input {2} -c 3 -i {3}
"""

def _select_annotations(input_directory,output_directory, annotation_header):
    header = open(annotation_header).readline().strip()
    logging.info('Generating combinations of {0} annotations'.format(len(header)))



def run_paintor(input_directory, annotation_header, auto_select_annotations = False,output_directory=None):
    """
        Function runs PAINTOR and selections the annotations for using Downstream.
    """    
    annotations = None 
    if auto_select_annotations:
        annotations = _select_annotations(input_directory, output_directory, annotation_header)
    if annotations is None:

    else:
    
    
