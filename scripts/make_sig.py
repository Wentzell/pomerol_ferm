#!/usr/bin/python

import h5py
import numpy as np
import sys
import os
import subprocess

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

# --- Constants
BETA=20.0

# --- Create G
gwimag = np.loadtxt("gw_imag00.dat")        # Load file into numpy array 
fdim_g = 2* len( set( gwimag[:,0] ) )
print "Counted " + str(fdim_g) + " fermionic frequencies for giw "
data_g = gwimag[:,1] + 1j * gwimag[:,2]

arr_g = 1j *np.zeros((fdim_g, 2, 2))
for i in range(fdim_g/2):
    arr_g[ fdim_g/2 + i,0,0] = data_g[i]
    arr_g[ fdim_g/2 + i,1,1] = data_g[i]

    arr_g[ fdim_g/2 - i - 1,0,0] = data_g[i].conjugate()
    arr_g[ fdim_g/2 - i - 1,1,1] = data_g[i].conjugate()

def G(w,i,j):
    return arr_g[fdim_g/2-fdim/2+w,i,j]

# --- Create Sig
gw0imag = np.loadtxt("gw0_imag00.dat")        # Load file into numpy array 
fdim_g = 2* len( set( gw0imag[:,0] ) )
print "Counted " + str(fdim_g) + " fermionic frequencies for giw0 "
data_g0 = gw0imag[:,1] + 1j * gw0imag[:,2]

arr_g0 = 1j *np.zeros((fdim_g, 2, 2))
for i in range(fdim_g/2):
    arr_g0[ fdim_g/2 + i,0,0] = data_g0[i]
    arr_g0[ fdim_g/2 + i,1,1] = data_g0[i]

    arr_g0[ fdim_g/2 - i - 1,0,0] = data_g0[i].conjugate()
    arr_g0[ fdim_g/2 - i - 1,1,1] = data_g0[i].conjugate()

arr_sig = 1j *np.zeros((fdim_g, 2, 2))

arr_sig[:,0,0] = 1. / arr_g0[:,0,0] - 1. / arr_g[:,0,0]
arr_sig[:,1,1] = 1. / arr_g0[:,1,1] - 1. / arr_g[:,1,1]

reG = arr_g.real
imG = arr_g.imag

reSig = arr_sig.real
imSig = arr_sig.imag

## Create and write to hdf5 file
print "Writing binary data to full.h5 using single precision and gzip compression"
f = h5py.File("full.h5", "w")

Ggrp = f.create_group("Giw")

dset = Ggrp.create_dataset("RE", (fdim_g,1,2,2), dtype='float32', compression="gzip", compression_opts=4) 
dset[:,0,:,:] = reG

dset = Ggrp.create_dataset("IM", (fdim_g,1,2,2), dtype='float32', compression="gzip", compression_opts=4) 
dset[:,0,:,:] = imG

Siggrp = f.create_group("Sig")

dset = Siggrp.create_dataset("RE", (fdim_g,1,2,2), dtype='float32', compression="gzip", compression_opts=4) 
dset[:,0,:,:] = reSig

dset = Siggrp.create_dataset("IM", (fdim_g,1,2,2), dtype='float32', compression="gzip", compression_opts=4) 
dset[:,0,:,:] = imSig

#print "Done!"
