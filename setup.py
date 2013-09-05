import os, sys
from glob import glob
#from distutils.core import setup,Extension

VERSION = "0.2.dev"
ISRELEASED = False
__author__ = "Gurpreet"
__version__ = VERSION

# metadata for setup()
metadata = {
    'version': VERSION,
    'author': __author__,
    'license': 'GPL version 2 or later',
    'install_requires': ['numpy' ],
    'platforms': ["Linux"],
    'zip_safe': False,
    'description': "Tools for data analysis",
    'long_description': """The code contains tools that were useful in my research.
    Most of the code in this collection is either copied directly or 
    is inspired by other codes, that I had to tinker with, to make it more suitable for my work
    You are free to do the same :-)
    """}





# setuptools needs to come before numpy.distutils to get install_requires
import setuptools 
import numpy
from distutils import sysconfig
from numpy.distutils.core import setup, Extension
from numpy.distutils.misc_util import Configuration

def configuration():
    config = Configuration('gp_grompy',package_parent="",top_path=None,
                           package_path='src/python')

    config.set_options(assume_default_configuration=True, delegate_options_to_subpackages=True,
                       quiet=False)
    # add the scipts, so they can be called from the command line
    #config.add_scripts([e for e in glob('scripts/*.py') if not e.endswith('__.py')])
    
    # add scripts as a subpackage (so they can be imported from other scripts)
    #config.add_subpackage('scripts',subpackage_path=None)

    # add gp_wham subpackage
    #config.add_subpackage('wham', subpackage_path='src/python/wham')


    return config

metadata['configuration'] = configuration
setup(**metadata)



