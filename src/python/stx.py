
##############################################################################3
# NOTE: Work in progress. The code is built on Martin hoefling's code 
# Visit git://github.com/GromPy/GromPy.git for further details.
# Author: Gurpreet Singh
# Date: 12-Feb-13
# version 0.1
###############################################################################
from ctypes import c_char, byref, POINTER, c_int, c_char_p
import ctypes
from gp_grompy.types import t_topology, t_atoms,t_trxstatus,output_env_t,\
                         t_trxframe,gmx_rmpbc_t,t_fileio
from gp_grompy import libgmx, rvec, matrix, c_real
import numpy as np
import math,sys,os
import time,datetime
# from grompy.tpxio import *
import logging
logger = logging.getLogger(__name__)
#logger = logging.getLogger('gmxwrapper')
#logger.setLevel(logging.DEBUG)


class Gmstx():
    def __init__(self):
        self.title = (c_char * 256)()
        self.top = t_topology()
        self.xp = POINTER(rvec)()
        self.vp = POINTER(rvec)()
        self.cePBC = c_int(-1)
        self.box = matrix()
        self.bMass = c_int()
    
    def read_tpr(self, filename):
        if not os.path.isfile(filename):
            raise SystemExit(filename,"not found")
        libgmx.init_top(byref(self.top))
        ret = libgmx.read_tps_conf(filename, byref(self.title), byref(self.top), \
        byref(self.cePBC), byref(self.xp), byref(self.vp), self.box, byref(self.bMass))
        if ret != 1:
            raise SystemExit("Error reading topology with wrapper call to read_tps_conf")
        logger.debug("%s File read", filename)
        self.natoms=self.top.atoms.nr
        self.atoms=self.top.atoms
        self.v=self.vp[0]
        self.fmt='tpr'
        
    def read_stx(self, filename):
        if not os.path.isfile(filename):
            raise SystemExit(filename,"not found")
        self.atoms = t_atoms()
        self.natoms = c_int()
        libgmx.get_stx_coordnum(filename, byref(self.natoms))    
        self.atoms.nr = c_int(self.natoms.value)
        self.xp = (rvec * self.natoms.value)()
        self.v = (rvec * self.natoms.value)()
        libgmx.init_t_atoms(byref(self.atoms), self.atoms.nr, c_int(0))
        libgmx.read_stx_conf(filename, self.title, byref(self.atoms), byref(self.xp), \
        byref(self.v), byref(self.cePBC), self.box)
        self.natoms=self.natoms.value
        logger.debug("%s succesfuly read", filename)
        self.fmt='nontpr'   
        
    def write_stx(self, filename,pbc=False):
        '''Writing to gro file is only implimented at the moment
           Also cannot write from tprfile 
        '''
        natoms=c_int(self.natoms)
        if hasattr(self,'fmt'):
            if self.fmt == 'tpr':
                raise SystemExit("writing from tpr not implimented")
        if pbc==True:
            gpbc=libgmx.gmx_rmpbc_init(byref(self.top.idef),self.cePBC,natoms,self.box)
            libgmx.gmx_rmpbc(gpbc,self.natoms,self.box,byref(self.xp))
            
        libgmx.write_sto_conf(filename, self.title, byref(self.atoms), byref(self.xp),
        byref(self.v), self.cePBC, self.box)    
        logger.debug("%s succesfuly written", filename)
        
    def copy(self, G):
        self.title = G.title
        self.atoms = G.atoms
        self.xp = G.xp
        self.v = G.v
        self.cePBC = G.cePBC
        self.box = G.box
        self.natoms=G.natoms
        self.fmt=G.fmt
        
    def x_to_array(self):
        ''' convert coordinates to numpy array
        '''
        natoms=self.natoms
        fr=np.zeros((natoms,3),dtype=np.float32)
        for i in xrange(natoms):
            for j in xrange(3):
                fr[i,j]=self.xp[i][j]
        
        return fr

    def array_to_x(self,arr):
        ''' convert coordinates to numpy array
        '''
        #self.x=arr
        pt=(rvec*self.natoms)
        self.xp=pt()
        for i in xrange(self.natoms):
            for j in xrange(3):
                self.xp[i][j]=arr[i,j]
        
        #self.xp=self.x.ctypes.data_as(pt)
