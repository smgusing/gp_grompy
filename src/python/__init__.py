################################################################
# NOTE: Work in progress. The code is built on Martin hoefling's code 
# Visit git://github.com/GromPy/GromPy.git for further details.

## The file is a modified copy of grompy.
################################################################
import os
import glob
import sys
from os import environ
import logging
import pkg_resources as pkgres
### Set up version information####
__version__ = pkgres.require("gp_grompy")[0].version


__all__ = ["stx","xtc","ndx","types"]

########################################
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
sh = logging.StreamHandler(sys.stdout)
#formatter = logging.Formatter(fmt='[%(name)s] %(message)s')
formatter = logging.Formatter(fmt='[%(asctime)s %(name)s] %(message)s',datefmt='%I:%M:%S')
sh.setFormatter(formatter)
logger.addHandler(sh)
logger.propagate = False
#########################################


from ctypes import c_float,\
                   cdll,\
                   c_double,\
                   c_int,\
                   pythonapi,\
                   py_object,\
                   c_void_p,POINTER

if environ.has_key("GMXLDLIB"):
    pth=environ["GMXLDLIB"]
else: 
    raise SystemExit("env variable GMXLDLIB not set")


logger.info("Loading gp_grompy with single precision library...")
c_real   = c_float
libmdname="%s/%s"%(pth,"libmd.so")
libgmxname="%s/%s"%(pth,"libgmx.so")
libmdrunname="%s/%s"%(pth,"mdrun.so")
#print libmdname
#libc     = cdll.LoadLibrary(libcname)
libmd    = cdll.LoadLibrary(libmdname)
libgmx   = cdll.LoadLibrary(libgmxname)
#libmdrun = cdll.LoadLibrary(libmdrunname)

# gromacs should be compiled with -lfftw3f!!!
#    libmdrun = cdll.LoadLibrary(environ["HOME"]+"/src/gromacs-4.0.5_TEST/src/kernel/libmdrun0/mdrun.so") # for debugging purposes


#FILE * to stderr, stdout
pythonapi.PyFile_AsFile.restype = c_void_p
stderr    = pythonapi.PyFile_AsFile(py_object(sys.stderr))
stdout    = pythonapi.PyFile_AsFile(py_object(sys.stdout))

rvec      = c_real*3
dvec      = c_double*3
ivec      = c_int*3
splinevec = c_real*3
matrix    = c_real*3*3
tensor    = c_real*3*3
########################################################################
from ndx import Gmndx
from xtc import Gmxtc
from stx import Gmstx
########################################################################









