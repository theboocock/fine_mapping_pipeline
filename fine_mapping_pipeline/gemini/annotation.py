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
import logging
import os
from gemini import GeminiQuery

class Annotation:
    """
        Represents one class of annotation for use in paintor.

        TODO extend these annotations to use general encode not just GEMINI
    """

    def __init__(self):
        self.annotation_rows = []
        self.unique_list = set()
    def add_data(self, key ,data):
        """
            Deal with the special cases of data from the gemini database
        """
        if ('encode_dnaseI_cell_list' == key):
            self.annotation_rows.append([(key, data)])
            if data is not None:
                data = data.split(';')
                for d in data:
                    self.unique_list.add(d)
        elif data != 'unknown' and data is not None:
            # Must be a cell line annotation so let's create that
            self.annotation_rows.append([(key, key + '_' + data)])
            self.unique_list.add(key + '_' +  data)
        else:
            self.annotation_rows.append([(key, data)])
    def get_annotation_matrix(self, row_num):
        """
            Returns a matrix for each annotation
        """
        matrix = []
        for u in self.unique_list:
            for r in self.annotation_rows[row_num]:
                key = r[0]
                value = r[1]
                if value is not None and value != 'unknown':
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
                        if value == u:
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
        annotation_matrix = []
        for i in range(self.snp_number):
            annotation_row = []
            for val in self.annotation_dict.values():
                annotation_row.extend(val.get_annotation_matrix(i))
            annotation_matrix.append(annotation_row)
        return annotation_matrix

    def get_header(self):
        """
            Return the header line for the annotation.
        """
        header = []
        for val in self.annotation_dict.values():
            header.extend(list(val.unique_list))
        return header

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
        return self.annotation_list.get_annotation_matrix()
    def get_snp_annotation_header(self):
        return self.annotation_list.get_header()


def generate_and_write_encode_annotations(loci, databases, output_directory):
    """
        Generate annotations for each locus indendently.

        @param databases gemini databases
        @param output_directory output dir
    """
    header_set = []
    snp_annotations = SNPAnnotations()
    snp_cutoffs =  []
    i = 0
    for database in databases:
        gemi_query = GeminiQuery(database)
        gemi_query.run('select encode_dnaseI_cell_list, encode_consensus_gm12878, encode_consensus_h1hesc, encode_consensus_helas3, encode_consensus_hepg2, encode_consensus_huvec, encode_consensus_k562 from variants')
        for row in gemi_query:
            snp_annotations.process_row(row)
            i += 1
        snp_cutoffs.append(i)
    header = snp_annotations.get_snp_annotation_header()
    matrix = snp_annotations.get_snp_annotation_matrix()
    i = 0
    j = 0
    for locus, snp_cutoff in zip(loci, snp_cutoffs):
        with open(os.path.join(output_directory, locus.rsid + '.annotations'), 'w') as out_annot:
            out_annot.write(str(snp_cutoff-i)+ ' ' + str(len(header)) + '\n')
            #out_annot.write(str(snp_cutoff-i)+ ' ' + "1" + '\n')
            for j in range(i, snp_cutoff):
                tmp_row = matrix[j]
                out_annot.write(' '.join([str(o) for o in  tmp_row]) + '\n')
            #    out_annot.write(str(tmp_row) + '\n')
        i = j + 1
    with open(os.path.join(output_directory, 'annotation.header'), 'w') as header_f:
        header_f.write(' '.join(header) + '\n')





        
