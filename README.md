# Fine mapping pipeline.

[![Build Status](https://travis-ci.org/smilefreak/fine_mapping_pipeline.svg?branch=master)](https://travis-ci.org/smilefreak/fine_mapping_pipeline)

Prepares summary statistics Z Score Files from IMPG and runs a finemapping analysis
using Paintor, Caviar, and BIMBAM (Caviarbf) 


## Installation and Setup.

Run the following commands to generate the files

### Prerequesites.

    - numpy
    - scipy
    - GEMINI framework - to install follow the instructions on the following github page https://github.com/arq5x/gemini
    
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

