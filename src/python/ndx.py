
from ctypes import c_char_p, c_int, pointer, POINTER, byref
from gp_grompy import libgmx
import logging
logger = logging.getLogger(__name__)

class Gmndx():
    
 def read_index(self, statfile, ngroups, isizep=None, indexp=None, grpnamesp=None):
    filename = c_char_p(statfile)
    if not isizep:
        self.isize = (c_int * ngroups)()
        isizep = pointer(self.isize)
    if not indexp:
        self.index = (POINTER(c_int) * ngroups)()
        indexp = pointer(self.index)
    if not grpnamesp:  
        self.grpnames = (c_char_p * ngroups)()
        grpnamesp = pointer(self.grpnames)
    
    libgmx.rd_index(filename, ngroups, isizep, indexp, grpnamesp)
  
    # return isize,index,grpnames
  
  
 def get_index(self, atoms, ngroups, filename=c_char_p(), isizep=None, indexp=None, grpnamesp=None):
    if type(filename) == str:
        filename = c_char_p(filename)

    if not isizep:
        self.isize = (c_int * ngroups)()
        isizep = pointer(self.isize)
    if not indexp:
        self.index = (POINTER(c_int) * ngroups)()
        indexp = pointer(self.index)
    if not grpnamesp:  
        self.grpnames = (c_char_p * ngroups)()
        grpnamesp = pointer(self.grpnames)

    libgmx.get_index(byref(atoms), filename, ngroups, isizep, indexp, grpnamesp)


