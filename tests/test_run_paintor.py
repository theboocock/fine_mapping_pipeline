import pytest

from fine_mapping_pipeline.finemap.paintor import run_paintor, _select_annotations, _do_lrt, _get_likelihood 
import logging 
logging.basicConfig(level=logging.INFO)

def test_run_paintor(tmpdir):
    input_directory = 'tests/paintor_data/'
    output_directory = tmpdir.mkdir('output')
    run_paintor(input_directory, annotation_header='tests/paintor_data/annotation.header',
                output_directory=output_directory)
    # TODO Check and read the data from the directory. 

def test_select_annotations(tmpdir):
    input_directory = 'tests/paintor_data/'
    annotations = _select_annotations(input_directory, causal_snp_number=3,
                                      annotation_header='tests/paintor_data/annotation.header')
    assert annotations ==  [1]

def test_do_lrf():
    p_value = _do_lrt(-10039, -10036)
    assert p_value == 0.014305878435429631 

def test_get_likelihood(tmpdir):
    input_directory = 'tests/paintor_data/'
    annot1 = _get_likelihood(input_directory, 1, 'TESTANNOTATION1', 3)
    assert annot1 == -279.196775
    annot2 = _get_likelihood(input_directory, 0, 'TESTANNOTATION2', 3)
    assert annot2 == -284.941121
