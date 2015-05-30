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

import argparse
from gemini import GeminiQuery

class Annotation:
    """
        Represents one class of annotation for use in paintor.
    """

    def __init__(self):
        self.annotation_rows = []
        self.unique_list = set()
    def add_data(self, key ,data):
        """

            Deal with the special cases of data from gemini

        """
        self.annotation_rows.append([(key,data)])
        if ('encode_dnaseI_cell_list' == key):
            if data is not None:
                data = data.split(';')
                for d in data:
                    self.unique_list.add(d)

        elif data is not None:
            self.unique_list.add(data)

    def get_annotation_matrix(self, row_num):
        """
            Returns a matrix for each annotation
        """
        matrix = []
        for u in self.unique_list:
            for r in self.annotation_rows[row_num]:
                key = r[0]
                value = r[1]
                if ( value is not None):
                    if(';' in value):
                        one_match = False
                        r_s = value.split(';')
                        for rs in r_s:
                            if rs == u:
                                matrix.append(1)
                                one_match = True
                                break
                        if not one_match:
                            matrix.append(0)

                    else:
                        if r == u:
                            matrix.append(1)
                        else:
                            matrix.append(0)
                else:
                    matrix.append(0)
        return matrix

class AnnotationMatrix:
    """
        Annotation Matrixx Object

        This class represents all the annotations that are used
    """

    def __init__(self):
        self.annotation_dict = {}
        self.snp_number =  0

    def add_annotation(self, key, data, snp_number):
        """
            Creates an annotation for use with paintor
        """
        try:
            self.annotation_dict[key].add_data(key, data)
        except KeyError:
            self.annotation_dict[key] = Annotation()
            self.annotation_dict[key].add_data(key, data)
        self.snp_number = snp_number

    def get_annotation_matrix(self):
        """
            Process each row individually

        """
        annotation_matrix = ''
        for i in range(self.snp_number):
            annotation_row = []
            for val in self.annotation_dict.values():
                annotation_row.extend(val[1].get_annotation_matrix(i))
            annotation_matrix += ' '.join([str(o) for o in annotation_row]) + '\n'
        return annotation_matrix

    def get_header(self):
        """
            Return the header line for the annotation.
        """
        header = []
        for val in self.annotation_dict.values()[0][0]:
            header.append(val)  
        return ' '.join(header)

class SNPAnnotations:
    """
        Represents annotations for all the SNPS.
    """
    def __init__(self):
        self.annotation_list = AnnotationMatrix()
        self.snp_number = 0

    def process_row(self, row):
        self.snp_number += 1
        for item in row:
            annotation_list = item
            break
        for annot in annotation_list:
            self.annotation_list.add_annotation(annot, row[annot], self.snp_number)

    def get_snp_annotation_matrix(self):
        self.annotation_list.get_annotation_matrix()
    def get_snp_annotation_header(self):
        self.annotation_list.get_header()


def generate_encode_annotations(databases):
    """
        Generate annotations for each locus indendently.
    """
    header_set = []
    for database in databases:
        gq.run('select encode_dnaseI_cell_list, encode_consensus_gm12878, encode_consensus_h1hesc, encode_consensus_helas3, encode_consensus_hepg2, encode_consensus_huvec, encode_consensus_k562 from variants')
        snp_annotations = SNPAnnotations()
        i = 0
        for row in gq:
            snp_annotations.process_row(row)
        snp_annotations.get_snp_annotation_matrix()
        header = snp_annotations.get_snp_annotation_header()
