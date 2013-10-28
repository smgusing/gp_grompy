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
    'author_email': 'togurpreet@gmail.com',
    'url': 'https://sites.google.com/site/togurpreet/Home',
    'install_requires': ['numpy' ],
    'platforms': ["Linux"],
    'zip_safe': False,
    'description': "python wrappers for gromacs",
    'long_description': """
    Code is based on grompy.
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
    return config

metadata['configuration'] = configuration
setup(**metadata)



