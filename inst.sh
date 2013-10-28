#! /bin/bash
source $HOME/.bashrc
lpy
##################################################################################
install_dir=/home/gurpreet/sft/python/env1
#rm -rvf build/ dist/ pmfcalculator.egg-info/
#rm -rvf /home/gurpreet/sft/python/env1/lib/python2.7/site-packages/pmfcalculator-0.1.dev-py2.7.egg
yes | pip uninstall gp_grompy
python setup.py install --prefix $install_dir
