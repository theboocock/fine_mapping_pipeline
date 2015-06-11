# Fine mapping pipeline.


Prepares summary statistics Z Score Files from IMPG and runs a finemapping analysis
using Paintor and Caviarbf.


Paintor. 
[https://github.com/gkichaev/PAINTOR_FineMapping/] 
Kichaev, Gleb, et al. "Integrating functional data to prioritize causal variants in statistical fine-mapping studies." (2014): e1004722.

Caviarbf.
[https://bitbucket.org/Wenan/caviarbf]
Chen, Wenan, et al. "Fine Mapping Causal Variants with an Approximate Bayesian Method Using Marginal Test Statistics." Genetics (2015): genetics-115.


## Installation and Setup.

Run the following commands to generate the files

### Prerequesites.

    - numpy
    - scipy
    - GEMINI framework [https://github.com/arq5x/gemini/]
    - PAINTOR
    - caviarbf
    - tabix [http://www.htslib.org/doc/tabix.html]
    - 
    
```

python setup.py install
bash install_paintor.sh

```

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

