#!/bin/sh
#
# description: get everything set up after logging on
#

# default python 3.11.4 on scinet's trillium cluster
# python -m venv .venv #make virtual env
source .venv/bin/activate #load virtual env
# pip install -r requirements.txt #install libraries
# install tsconvert, not on pypi
# git clone https://github.com/tskit-dev/tsconvert.git
# cd tsconvert
# pip install .
# cd -

# get relate
# git clone https://github.com/MyersGroup/relate.git
# cd relate/build
# cmake ..
# make
# cd -
# note: cmake version 3.27.7, gcc version 12.3.1 20230526

export XDG_CACHE_HOME=$SCRATCH #so snakemake writes temp files where it is allowed

# commented out lines 62-65 of .venv/lib/python3.11/site-packages/snakemake_executor_plugin_slurm/submit_string.py to prevent mem requests to slurm (not allowed on trillium)
