from setuptools import setup, find_packages

import os 

def locate_packages():
    packages = ['fine_mapping_pipeline']
    for (dirpath, dirnames, _) in os.walk(packages[0]):
        for dirname in dirnames:
            package = os.path.join(dirpath, dirname).replace(os.sep, ".")
            packages.append(package)
    return packages

setup(
    name="fine_mapping_pipeline",
    version="1.0",
    packages=locate_packages(),
    author="James Boocock",
    author_email="james.boocock@otago.ac.nz",
    description="PAINTOR, Caviar and BIMBAM pipeline for fine mapping",
    license="Mit",
    zip_safe=False,
     entry_points={
        'console_scripts': [
            'fine_mapping_pipeline= fine_mapping_pipeline.pipeline:main',
        ]
        },
    url="github.com/smilefreak/fine_mapping_pipelin",
    use_2to3=True,
)
