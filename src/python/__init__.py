################################################################
## The file is a modified copy of grompy.
## See Grompy for original author information
## Modified By, Gurpreet Singh
## version 0.10
## 
################################################################
import os
import glob
import sys
from os import environ


__all__ = ["stx","xtc","ndx","types"]

import logging
logger = logging.getLogger(__name__)
#logger.setLevel(logging.DEBUG)
#ch = logging.StreamHandler(sys.stdout)
#formatter = logging.Formatter(fmt='%(filename)s:%(funcName)s:%(message)s')
#ch.setFormatter(formatter)
#logger.addHandler(ch)
#logger.propagate = False


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
#from rmsfit import Rmsfit
#import types
########################################################################









