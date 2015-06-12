# Fine mapping pipeline.


Prepares summary statistics Z Score Files from IMPG and runs a finemapping analysis
using Paintor and Caviarbf.

## Citations

If you use this tool please ensure you cite the primary papers for the methods, these are listed below. 

###Paintor

https://github.com/gkichaev/PAINTOR_FineMapping/
 
Kichaev, Gleb, et al. "Integrating functional data to prioritize causal variants in statistical fine-mapping studies." (2014): e1004722.

###Caviarbf

https://bitbucket.org/Wenan/caviarbf

Chen, Wenan, et al. "Fine Mapping Causal Variants with an Approximate Bayesian Method Using Marginal Test Statistics." Genetics (2015): genetics-115.


## Installation and Setup.

Run the following commands to generate the files

### Prerequesites.

    - numpy
    - scipy
    - GEMINI framework https://github.com/arq5x/gemini/
    - PAINTOR https://github.com/gkichaev/PAINTOR_FineMapping/
    - caviarbf https://bitbucket.org/Wenan/caviarbf
    - tabix http://www.htslib.org/doc/tabix.html

#### IMPORTANT

Before installation it is important you change the file fine_mapping_pipeline/config/__init__.py to either a local installation
or the 1000 genomes FTP server. The {0} is replaced internall for each chromosome. The following string is an example for the EBI 1000 genomes FTP.

"ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{0}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

```

python setup.py install
bash install_paintor.sh

```

## Running the pipeline. 

Basically the pipeline takes a list of either rsids, or chromosome"\t"position, and a folder containing 
ZScores statistics that have been output using the ImpG program, or any program as long as it conforms to the output we expect. 

https://github.com/huwenboshi/ImpG

The pipeline looks at the 1st two columns for chromosome and position, and the 5th column for the ZScore. 

Example ZScore input files are found in the tests/test_data/zscores/ folder.
Example RSID file for the Zscores is found in the tests/test_data/rsids.txt 

To get familiar with the pipeline run the following commands on the test dataset.

First you need to prepare you ImpG data for running the pipeline.

```
    fine_mapping_pipeline prepare -s tests/test_data/rsids.txt -z tests/test_data/zscores/ -f 25000 -b hg19 -o tests/test_output
```

Then follow up with either paintor and/or caviarbf with you parameters of interest, more customisation to be added soon.

Run these commands to get a feel.

``` 
    fine_mapping_pipeline paintor -d tests/test_output -o tests/test_output/paintor -c 1
    fine_mapping_pipeline caviarbf -d tests/test_output -o tests/test_output/caviarbf -c 1
```

This will place the analysis results from paintor and caviarbf in the tests/test_output/paintor and tests/test_output/caviarbf directories.

In the future, I plan to add plotting tools of the results data. 

Also, currently PAINTOR uses a bonferonni correction for the number of annotations used, to restrict the number of annotations to say the 
genome segmentations use the argument --no-dhs to the paintor sub-command of the fine-mapping pipeline.

## Testing

Tox is used to perform Unit testing of the Pipeline, make sure it is installed
and accessible on your command-line.

```

tox

```

or alteratively you can use py.test directly.

```

py.test

```

