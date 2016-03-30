#!/usr/bin/python

import h5py
import numpy as np
import sys
import os
import subprocess

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

run('./merge.sh')

# --- Constants
BETA=10.0

print "Filling arrays ..."

#data1111 = np.loadtxt("1111.dat")        # Load file into numpy array 
data1010 = np.loadtxt("1010.dat")        # Load file into numpy array 
#data1001 = np.loadtxt("1001.dat")        # Load file into numpy array 

# Count fermionic frequencies
fdim = len( set( data1010[:,0] ) )
print "Counted " + str(fdim) + " fermionic frequencies for vertex "

# Prepare data for output
arr1010 = np.reshape(data1010[:,3], (fdim, fdim, fdim)) + 1j * np.reshape(data1010[:,4], (fdim, fdim, fdim))

# --- Create G
gwimag = np.loadtxt("gw_imag00.dat")        # Load file into numpy array 
fdim_g = 2* len( set( gwimag[:,0] ) )
print "Counted " + str(fdim_g) + " fermionic frequencies for giw "
data_g = gwimag[:,1] + 1j * gwimag[:,2]

arr_g = 1j *np.zeros((fdim_g))
for i in range(fdim_g/2):
    arr_g[ fdim_g/2 + i ] = data_g[i]
    arr_g[ fdim_g/2 - i - 1 ] = data_g[i].conjugate()

def G(w):
    return arr_g[fdim_g/2-fdim/2+w]

# --- Create Sig
gw0imag = np.loadtxt("gw0_imag00.dat")        # Load file into numpy array 
fdim_g = 2* len( set( gw0imag[:,0] ) )
print "Counted " + str(fdim_g) + " fermionic frequencies for giw0 "
data_g0 = gw0imag[:,1] + 1j * gw0imag[:,2]

arr_g0 = 1j *np.zeros((fdim_g))
for i in range(fdim_g/2):
    arr_g0[ fdim_g/2 + i ] = data_g0[i]
    arr_g0[ fdim_g/2 - i - 1 ] = data_g0[i].conjugate()

arr_sig = 1. / arr_g0[:] - 1. / arr_g[:]

# -- Vert
print("Calculating vert ...")

# Subtract free propagation 

for w1 in range(fdim):
    for w2 in range(fdim):
        arr1010[w1,w2,w1] -= BETA * G(w1) * G(w2) # subtract non-crossing diagram

## Assume SU2 symmetry for now
for w1 in range(fdim):
    for w2 in range(fdim):
        for w1p in range(fdim):
            w2p = w1+w2-w1p
            arr1010[w1,w2,w1p] /= G(w1)*G(w2)
            arr1010[w1,w2,w1p] /= G(w1p)*G(w2p)

reVert = arr1010.real
imVert = arr1010.imag

reG = arr_g.real
imG = arr_g.imag

reSig = arr_sig.real
imSig = arr_sig.imag

## Create and write to hdf5 file
print "Writing binary data to full.h5 using single precision and gzip compression"
f = h5py.File("full.h5", "w")

Vertgrp = f.create_group("Vert")

dset = Vertgrp.create_dataset("RE", (fdim, fdim, fdim,1,1,1,1,1,1,1), dtype='float32', compression="gzip", compression_opts=4) # compression_opts from 0 - 9(strong compression) # or dtype='float64' for double precision
dset[:,:,:,0,0,0,0,0,0,0] = reVert

dset = Vertgrp.create_dataset("IM", (fdim, fdim, fdim,1,1,1,1,1,1,1), dtype='float32', compression="gzip", compression_opts=4) # compression_opts from 0 - 9(strong compression) # or dtype='float64' for double precision
dset[:,:,:,0,0,0,0,0,0,0] = imVert

Ggrp = f.create_group("Giw")

dset = Ggrp.create_dataset("RE", (fdim_g,1,1,1), dtype='float32', compression="gzip", compression_opts=4) 
dset[:,0,0,0] = reG

dset = Ggrp.create_dataset("IM", (fdim_g,1,1,1), dtype='float32', compression="gzip", compression_opts=4) 
dset[:,0,0,0] = imG

Siggrp = f.create_group("Sig")

dset = Siggrp.create_dataset("RE", (fdim_g,1,1,1), dtype='float32', compression="gzip", compression_opts=4) 
dset[:,0,0,0] = reSig

dset = Siggrp.create_dataset("IM", (fdim_g,1,1,1), dtype='float32', compression="gzip", compression_opts=4) 
dset[:,0,0,0] = imSig

#print "Done!"
