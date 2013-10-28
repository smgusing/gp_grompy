#! /usr/bin/env python
# -*- coding: utf-8 -*-

##############################################################################3

# NOTE: Work in progress. The code is built on Martin hoefling's code 
# Visit git://github.com/GromPy/GromPy.git for further details.
# Author: Gurpreet Singh
# Date: 12-Feb-13
# version 0.1
###############################################################################
from ctypes import c_char, byref, POINTER, c_int, c_char_p
import ctypes
from gp_grompy.types import t_topology, t_atoms, t_trxstatus, output_env_t, \
                         t_trxframe, gmx_rmpbc_t, t_fileio, GMXctypesError
from gp_grompy import libgmx, rvec, matrix, c_real
import numpy as np
import math, sys, os
import time, datetime
# from grompy.tpxio import *
import logging
ch = logging.StreamHandler()
formatter = logging.Formatter(fmt='[%(levelname)s:%(filename)s] %(message)s')
ch.setFormatter(formatter)
logger=logging.getLogger()
#logger.setLevel(logging.INFO)
logger.setLevel(logging.DEBUG)
logger.addHandler(ch)



class Gmxtc():
    def __init__(self,logger=None):
        self.natoms = c_int()
        self.step = c_int()
        self.time = c_real()
        self.prec = c_real()
        self.bOK = c_int()
        self.box = matrix()
        self.xp = POINTER(rvec)()
    
    def open_xtc(self, filename, mode):
	if mode == 'r':
            if not os.path.isfile(filename):
	        raise SystemExit(filename, "not found")
        fileptr = c_char_p(filename)
        modeptr = c_char_p(mode)
        libgmx.open_xtc.restype = POINTER(t_fileio)
        self.xtcfh = libgmx.open_xtc(fileptr, modeptr)
        logger.debug("%s File opened in %s mode", filename, mode)
        
    def read_first_xtc(self):
        ret = libgmx.read_first_xtc(self.xtcfh, byref(self.natoms), byref(self.step), \
              byref(self.time), self.box, byref(self.xp), byref(self.prec), \
              byref(self.bOK))
        if ret != 1:
                raise GMXctypesError, "read_first_xtc did not return 1"
        logger.debug("First frame read successfully")
        return ret       
        
    def read_next_xtc(self):
        self.step = c_int(-1)
        self.time = c_real(-1)
        self.prec = c_real(-1)
        self.bOK = c_int()
        ret = libgmx.read_next_xtc(self.xtcfh, self.natoms, byref(self.step), \
              byref(self.time), self.box, self.xp, byref(self.prec), byref(self.bOK))
        
        return ret
    
    def read_timeframe(self, filename, time=None):
        if not os.path.isfile(filename):
            raise SystemExit(filename, "not found")

        self.xtcfh = libgmx.gmx_fio_open(c_char_p(filename), c_char_p("r"))
        
        logger.debug("%s File opened", filename)
        
        self.read_first_xtc()
        ### bug in reading frame using seektime when it is the first frame##
        if time == 0:
            self.close_xtc()
            return
            
        if (time is not None) :
            if (self.time.value != 0):
                logger.warn("Initial time is not zero %s. Will add to time"%round(self.time.value,3))
                time += round(self.time.value,3)
        else:
            raise SystemExit("provide timeframe to extract")
        
        # int xtc_seek_time(t_fileio *fio, real time, int natoms);
        ret = libgmx.xtc_seek_time(self.xtcfh, c_real(time), self.natoms)
            
        if ret != 0:
                raise SystemExit("Cannot find step %f time %f\n", stepno, time)
            
        ret = self.read_next_xtc()
        
        if round(self.time.value,3) != time:
            mes="Time loaded {0} does not match time requested {1}".format(self.time.value,time)
            raise SystemExit(mes)
        
#         if ret != 1:
#                 raise SystemExit("Cannot read step %f, time %f\n", stepno, time)
        logger.debug("Time %f loaded successfuly:Step %d", self.time.value,
                      self.step.value)
        
        self.close_xtc()
        
    def read_stepframe(self, filename, stepno=None):
        
        if not os.path.isfile(filename):
            raise SystemExit(filename, "not found")

        self.xtcfh = libgmx.gmx_fio_open(c_char_p(filename), c_char_p("r"))
        
        logger.debug("%s File opened", filename)
        
        self.read_first_xtc()
            
        if (stepno is not None):
            if (self.step.value != 0):
                logger.warn("Initial step is not zero %s. Will add to step"%self.time.value)
                stepno += self.step.value
            ret = libgmx.xtc_seek_frame(self.xtcfh, c_int(stepno), self.natoms)
            
        if ret != 0:
                raise SystemExit("Cannot find step %f time %f\n", stepno, time)
            
        ret = self.read_next_xtc()
        
        if ret != 1:
                raise SystemExit("Cannot read step %f, time %f\n", stepno, time)
        logger.debug("Time %f loaded successfuly:Step %d", self.time.value,
                      self.step.value)
        self.close_xtc()
        
    def load_traj(self, xtcfile, skip=1, bPBC=0, nframes=1,
                  isize=None, index=None, tprfile=None):
        """
        """
        stime = time.time()
        # gppbc = gmx_rmpbc_t()
        buffer_from_memory = ctypes.pythonapi.PyBuffer_FromMemory
        buffer_from_memory.restype = ctypes.py_object
        if skip < 1: raise SystemExit("skip cannot be less than 1")
        if bPBC:
            tpr = Gmstx()
            tpr.read_tprfile(tprfile)
            gpbc = libgmx.gmx_rmpbc_init(byref(self.top.idef), self.cePBC, self.natoms, self.box)
            bPBC = c_int(bPBC)

        self.open_xtc(xtcfile, 'r')
        ret = self.read_first_xtc()

        if isize == None:
            natoms = self.natoms.value  # read full traj
            index = np.arange(natoms)
        else:
            natoms = isize

        coords = []
        itr = 0
        vflag = False
        fr0time = round(self.time.value,3)
        logger.debug("Timestamp on first frame %s",fr0time)

        while ret == 1:
            if bPBC == 1:
                libgmx.gmx_rmpbc(gpbc, self.natoms, self.box, self.xp)
            ##http://bugs.python.org/issue10746
            ###More elegant but buggy
            ### Triggers runtime warning    
            #fr = np.ctypeslib.as_array((rvec * natoms).from_address(ctypes.addressof(self.xp.contents)))
            # 32 bit float =4 bytes
            buf = buffer_from_memory(self.xp, 4 * natoms * 3)
            fr=np.ndarray((natoms, 3),dtype=np.float32, order='C',
                         buffer=buf)
            fr = fr[index, :]
            ## python does not own the memory so make a copy
            frcopy = np.copy(fr)
            coords.append(frcopy)
            itr = itr + 1
            if itr % 500 == 0 : logger.debug("%d frames loaded from %s", itr, xtcfile)
            # if we know the delta time between frames then 
            # move to next frame
            if ((skip > 1) and (vflag == True)):
                seektime = c_real(self.time.value + deltatime * skip)
                ret = libgmx.xtc_seek_time(self.xtcfh, seektime, self.natoms)

            # Read next frame
            if ret > -1:
                ret = self.read_next_xtc()

            # if we are doing this first time    
            if ((skip > 1) and (vflag == False)):
                vflag = True
                #seektime = c_real(self.time.value + initvalue * (skip - 1))
                #ret = libgmx.xtc_seek_time(self.xtcfh, seektime, self.natoms)
                #if ret > -1:    
                fr1time = round(self.time.value,3)
                deltatime=fr1time-fr0time
                logger.debug("Using %s as timestep between frames",deltatime)

                seektime = c_real(fr0time + deltatime * skip)
                ret = libgmx.xtc_seek_time(self.xtcfh, seektime, self.natoms)
                if ret > -1:
                    ret = self.read_next_xtc()

        coordout = np.vstack(coords).reshape(itr, natoms, 3)
        runtime = datetime.timedelta(seconds=(time.time() - stime))
        logger.info("%d frames from %s loaded", itr, xtcfile)
        logger.debug("%s time required loading %d frames from %s", runtime, itr, xtcfile)
        return coordout
       
    
    def load_all_frames(self, xtcfile, skip=1, bPBC=0, nframes=1,
                  isize=None, index=None, tprfile='xa.tpr'):
        """
        """
        stime = time.time()
        # gppbc = gmx_rmpbc_t()
        buffer_from_memory = ctypes.pythonapi.PyBuffer_FromMemory
        buffer_from_memory.restype = ctypes.py_object
        boxes=[]
        times=[]
        frnos=[]
        
        if skip < 1: raise SystemExit("skip cannot be less than 1")
        if bPBC:
            tpr = Gmstx()
            tpr.read_tprfile(tprfile)
            gpbc = libgmx.gmx_rmpbc_init(byref(self.top.idef), self.cePBC, self.natoms, self.box)
            bPBC = c_int(bPBC)
        self.open_xtc(xtcfile, 'r')
        ret = self.read_first_xtc()
        prec=self.prec
        if isize == None:
            natoms = self.natoms.value  # read full traj
            index = np.arange(natoms)
        else: 
            natoms = isize
        coords = []    
        itr = 0
        vflag = False
        while ret == 1:
            if bPBC == 1:
                libgmx.gmx_rmpbc(gpbc, self.natoms, self.box, self.xp)
            ##http://bugs.python.org/issue10746
            ###More elegant but buggy
            ### Triggers runtime warning    
            #fr = np.ctypeslib.as_array((rvec * natoms).from_address(ctypes.addressof(self.xp.contents)))
            # 32 bit float =4 bytes
            buf = buffer_from_memory(self.xp, 4 * natoms * 3)
            fr=np.ndarray((natoms, 3),dtype=np.float32, order='C',
                         buffer=buf)
            fr = fr[index, :]
            ## python does not own the memory so make a copy
            frcopy = np.copy(fr)
            coords.append(frcopy)
            boxes.append(self.box)
            times.append(self.time)
            itr = itr + 1
            if itr % 500 == 0 : logger.debug("%d frames loaded from %s", itr, xtcfile)
            if ((skip > 1) and (vflag == True)):
                seektime = c_real(self.time.value + initvalue * skip)
                ret = libgmx.xtc_seek_time(self.xtcfh, seektime, self.natoms)
            if ret > -1:    
                ret = self.read_next_xtc()
            if ((skip > 1) and (vflag == False)):
                vflag = True
                initvalue = self.time.value
                seektime = c_real(self.time.value + initvalue * (skip - 1))
                ret = libgmx.xtc_seek_time(self.xtcfh, seektime, self.natoms)
                if ret > -1:    
                    ret = self.read_next_xtc()
        coordout = np.vstack(coords).reshape(itr, natoms, 3)
        runtime = datetime.timedelta(seconds=(time.time() - stime))
        logger.info("%s time required loading %d frames from %s\n", runtime, itr, xtcfile)
        return coordout,boxes,times,prec

    def get_trajlength(self, xtcfile):
        """quick way to get number of frames in trajectory
            will only work if trajectories are continuous,frames have constant time diff
            and precision is 1000
        """
        
        stime = time.time()
        self.open_xtc(xtcfile, 'r')
        ret = self.read_first_xtc()
        natoms = self.natoms.value
        fr0time = self.time.value
        itr = 0
        vflag = False
        while ret == 1:
            itr = itr + 1
            ret = self.read_next_xtc()
            if itr % 500 == 0 : logger.debug("%d frames read from %s", itr, xtcfile)
                    
        runtime = datetime.timedelta(seconds=(time.time() - stime))
        logger.info("%s time required counting %d frames from %s\n", runtime, itr, xtcfile)
        return itr 


 
    def close_xtc(self):
        libgmx.close_xtc(self.xtcfh)

    def write_xtc(self):
        ret = libgmx.write_xtc(self.xtcfh, self.natoms, self.step, self.time, \
                                self.box, self.xp, self.prec)
        return ret
        
    def write_array_as_traj(self, filename, traj, boxs, times, prec):
        self.open_xtc(filename, 'w')
        logger.info("Opening %s for writing", filename)
        nframes, natoms = traj.shape[0], traj.shape[1]
        step = 0
        rvec_p = POINTER(rvec)
        # you are slave of two masters. Make sure not to release the memory
        # in python until you are done with it in C.
        for i in range(nframes):
            time = times[i]
            box = boxs[i]
            x = traj[i]
            xp = x.ctypes.data_as(rvec_p)
            step = i
            ret = libgmx.write_xtc(self.xtcfh, natoms, step, time, \
                                box, xp, prec)
             
        return ret

    def copy(self, G, index=None):
        self.step = G.step
        self.time = G.time
        self.box = G.box
        self.prec = G.prec
        
        if index != None:
            natoms = index.size
            self.natoms = c_int(natoms)
            rvec_p = POINTER(rvec)
            # you are slave of two masters. Make sure not to release the memory
            # in python until you are done with it in C.
            self.x = np.zeros((natoms, 3), dtype=np.float32, order='C')
            self.xp = self.x.ctypes.data_as(rvec_p)
            for i in xrange(natoms):
                for j in xrange(3):
                    self.xp[i][j] = G.xp[index[i]][j]
                    
        else:
            self.xp = G.xp
            self.natoms = G.natoms
            
    def x_to_array_slow(self):
        ''' convert coordinates to numpy array
        '''
        natoms = self.natoms.value
        fr = np.zeros((natoms, 3), dtype=np.float32)
        for i in xrange(natoms):
            for j in xrange(3):
                fr[i, j] = self.xp[i][j]
        
        return fr
    def x_to_array(self):
        ''' convert coordinates to numpy array
        '''
        buffer_from_memory = ctypes.pythonapi.PyBuffer_FromMemory
        buffer_from_memory.restype = ctypes.py_object
        #rvec_p = POINTER(rvec)
        natoms = self.natoms.value
        #fr = np.ctypeslib.as_array((rvec * natoms).from_address(ctypes.addressof(self.xp.contents)))
        buf = buffer_from_memory(self.xp, 4 * natoms * 3)
        fr=np.ndarray((natoms, 3),dtype=np.float32, order='C',
                     buffer=buf)
        ## python does not own the memory so make a copy
        frcopy = np.copy(fr)
        return frcopy

            
    def read_first_x(self, filename):
        print "not working yet"
        sys.exit(1)
        oenv = output_env_t()
        status = POINTER(t_trxstatus)()
        fileptr = c_char_p(filename)
        self.natoms = libgmx.read_first_x(oenv, byref(status), fileptr, byref(self.time), byref(self.xp), self.box)
        print self.natoms
