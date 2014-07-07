===============================================
gp_grompy
===============================================

:Author: Gurpreet Singh
:Contact: togurpreet@gmail.com
:License: Read LICENSE.txt 

-----------------------------------------------
Introduction
----------------------------------------------- 
gp_grompy is a python module that can be used to read and write
Gromacs files. The trajectories can be loaded into numpy arrays. The package draws its
inspiration from GromPy (https://github.com/GromPy/GromPy).
Likewise, it uses ctypes to interface with Gromacs libraries (libmd.so and libgmx.so)

The package is not as elaborate as GromPy as only very few Gromacs functions are wrapped,
mainly those that can read/write gromacs files and perform few operations such as pbc removal.

-----------------------------------------------
Installation
-----------------------------------------------

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Prerequisites
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
- Numpy 
- Gromacs (version 4.5.*) developed using 4.5.5
- gcc compiler
    + gcc-4.7.3 and above. 

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Quick install
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
python setup.py install


-------------------------------------------------
Usage
-------------------------------------------------
The environment variable ``GMXLDLIB`` should point to the directory
containing libmd.so and libgmx.so. This can be acheived by sourcing GMXRC from Gromacs
installation.

The directory ``test`` contains ``test.py``. The script shows the basic usage.

The package parallelclusterer uses gp_grompy to read xtc, tpr and ndx files.
Have a look at its code at https://github.com/smgusing/parallelclusterer for more usage examples 
   
------------------------------------------------
TO DO
------------------------------------------------
- The code is quite rudimentary. The Project needs structuring.
- Use cython instead of ctypes for interfacing with Gromacs. 
 









