from setuptools import setup, find_packages

import os 

def locate_packages():
    packages = ['paintor_pipeline']
    for (dirpath, dirnames, _) in os.walk(packages[0]):
        for dirname in dirnames:
            package = os.path.join(dirpath, dirname).replace(os.sep, ".")
            packages.append(package)
    return packages

setup(
    name="paintor_pipeline",
    version="1.0",
    packages=locate_packages(),
    author="James Boocock",
    author_email="james.boocock@otago.ac.nz",
    description="Use PAINTOR for any IMPG dataset",
    license="MIT",
    zip_safe=False,
     entry_points={
        'console_scripts': [
            'paintor_pipeline = paintor_pipeline.pipeline:main',
        ]
        },
    url="github.com/smilefreak/paintor_pipeline",
    use_2to3=True,
)
