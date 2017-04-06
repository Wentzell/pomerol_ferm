#!/usr/bin/python

import h5py
import numpy as np
import sys
import os
import subprocess
import optparse
import numpy.linalg as la

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

# --- Files for full and bare calculation
f = h5py.File("full.h5", "r")
f0 = h5py.File("bare.h5", "r") # no G2 needed

# --- Read data
parVals = f["/Params"].attrs.values()

UINT =  parVals[0] # follows order in output.cpp
BETA =  parVals[1]
Giw = np.array(f["/Giw/RE"]) + 1j*np.array(f["/Giw/IM"])
G2 = np.array(f["/G2/RE"]) + 1j*np.array(f["/G2/IM"])
Giw0 = np.array(f0["/Giw/RE"]) + 1j*np.array(f0["/Giw/IM"])

# --- Count fermionic frequencies
fdimG = Giw.shape[0]
fdimG2 = G2.shape[0]
qncount = Giw.shape[3]

shift_G = fdimG/2 - fdimG2/2

# --- Create Sig
Giw_inv = np.reshape( map( la.inv, np.squeeze(Giw) ), (fdimG,1,qncount,qncount) )
Giw0_inv = np.reshape( map( la.inv, np.squeeze(Giw0) ), (fdimG,1,qncount,qncount) )
Sig = Giw0_inv - Giw_inv

# -- Vert

print("Calculating vert ...")

G2c=G2

# Subtract free propagation 
for w1 in range(fdimG2):
    for w2 in range(fdimG2):
        for i in range(qncount):
            for j in range(qncount):
                for k in range(qncount):
                    for l in range(qncount): 

                        # subtract non-crossing diagram
                        G2c[w1,w2,w1,0,0,0,i,j,k,l] -= BETA * Giw[w1 + shift_G,0,i,k] * Giw[w2 + shift_G,0,j,l]

vert=1j*np.zeros(G2c.shape)

# Amputate legs
for w1 in range(fdimG2):
    print " .. Amputing for w1 " + str(w1)
    for w2 in range(fdimG2):
        for w1p in range(fdimG2):

            w2p = w1+w2-w1p

            # Loop over external indeces
            for i in range(qncount):
                for j in range(qncount):
                    for k in range(qncount):
                        for l in range(qncount):
                            
                            # Loop over internal indeces
                            for i1 in range(qncount):
                                for j1 in range(qncount):
                                    for k1 in range(qncount):
                                        for l1 in range(qncount):
                                             vert[w1,w2,w1p,0,0,0,i,j,k,l] += Giw_inv[w1 + shift_G,0,i,i1]*Giw_inv[w2 + shift_G,0,j,j1] * G2c[w1,w2,w1p,0,0,0,i1,j1,k1,l1] * Giw_inv[w1p + shift_G,0,k1,k]*Giw_inv[w2p + shift_G,0,l1,l]

reG2 = G2.real
imG2 = G2.imag

reVert = vert.real
imVert = vert.imag

reG = Giw.real
imG = Giw.imag

reSig = Sig.real
imSig = Sig.imag

# Create and write to hdf5 file
print "Writing binary data to final.h5 using single precision and gzip compression"
fout = h5py.File("final.h5", "w")

Pargrp = fout.create_group("Params")
Pargrp.attrs.create("UINT", UINT, dtype='float32')
Pargrp.attrs.create("BETA", BETA, dtype='float32')

G2grp = fout.create_group("G2")
G2grp.create_dataset("RE", data=reG2, dtype='float32', compression="gzip", compression_opts=4) # compression_opts from 0 - 9(strong compression) # or dtype='float64' for double precision
G2grp.create_dataset("IM", data=imG2, dtype='float32', compression="gzip", compression_opts=4) 
G2grp.create_dataset("fgrid", data=np.array(f["/G2/fgrid"]), dtype='float32') 

Vertgrp = fout.create_group("Vert")
Vertgrp.create_dataset("RE", data=reVert, dtype='float32', compression="gzip", compression_opts=4)
Vertgrp.create_dataset("IM", data=imVert, dtype='float32', compression="gzip", compression_opts=4)
Vertgrp.create_dataset("fgrid", data=np.array(f["/G2/fgrid"]), dtype='float32') 

Ggrp = fout.create_group("Giw")
Ggrp.create_dataset("RE", data=reG, dtype='float32') 
Ggrp.create_dataset("IM", data=imG, dtype='float32') 
Ggrp.create_dataset("fgrid", data=np.array(f["/Giw/fgrid"]), dtype='float32') 

Siggrp = fout.create_group("Sig")
Siggrp.create_dataset("RE", data=reSig, dtype='float32') 
Siggrp.create_dataset("IM", data=imSig, dtype='float32') 
Siggrp.create_dataset("fgrid", data=np.array(f["/Giw/fgrid"]), dtype='float32') 

print "Done!"
