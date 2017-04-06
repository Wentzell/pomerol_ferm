#!/usr/bin/python

#--------------------------------------IMPORTS ------------------------------------------

import h5py
import matplotlib.pyplot as pl
import numpy as np
import os
import sys
import subprocess
import math

#--------------------------------------SETTINGS ------------------------------------------
#if -1, change overall vertex sign (switch definition)
vert_mul = 1

shift=0

#--------------------------------------MANAGING FILES ------------------------------------------

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

fname = "full.h5"

if len(sys.argv) > 1:
    fname = str(sys.argv[1])

fname = fname.rstrip('\n') # strip newline of fname
f = h5py.File(fname, "r")

os.system('mkdir -p plots')

#--------------------------------------GENERAL PLOT SETTINGS------------------------------------------

pl.rc('xtick', labelsize=9) 
pl.rc('ytick', labelsize=9) 
#pl.rc('text', usetex=True)
#pl.rc('text.latex', preamble='\usepackage{amsmath}')

RE = r"$\operatorname{Re}"
IM = r"$\operatorname{Im}"

#--------------------------------------SELF ENERGY PLOTTING ------------------------------------------

print("Plotting self-energy ...")


#--- Read
resig = f["/Sig/RE"]
imsig = f["/Sig/IM"]

revert = vert_mul * np.array(f["/Vert/RE"])
imvert = vert_mul * np.array(f["/Vert/IM"])

fdim = resig.shape[0]

siggrid = f["/Sig/fgrid"]
vertgrid = f["/Vert/fgrid"]

#--- Helper functions
def plotSig( use_pl, arr, string ):
    pl.plot( siggrid, arr[:,0,0,0], 'bx', ms=5, mew=0.5)
    #pl.xlim([0,30])
    #pl.ylim([-0.08,0])
    use_pl.set_title(string)
    return


#--- Plot physical
plotSig( pl.subplot(2,2,1), resig, RE + "\Sigma(i\omega)$" ) 
plotSig( pl.subplot(2,2,2), imsig, IM + "\Sigma(i\omega)$" ) 

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--------------------------------------G(iw) PLOTTING ------------------------------------------

print("Plotting Green function ...")


#--- Read
reGiw = f["/Giw/RE"]
imGiw = f["/Giw/IM"]
fdim = reGiw.shape[0]

Giwgrid = np.arange(fdim)


#--- Helper functions
def plotGiw( use_pl, arr, string ):
    pl.plot( Giwgrid, arr[:,0,0,0], 'bx', ms=3, mew=0.2)
    #pl.xlim([1.2*min(vertgrid),1.2*max(vertgrid)])
    use_pl.set_title(string)
    return


#--- Plot physical
plotGiw( pl.subplot(2,2,1), reGiw, RE + "G(i\omega)$" ) 
pl.xlabel(r"$\omega_n$")
plotGiw( pl.subplot(2,2,2), imGiw, IM + "G(i\omega)$" ) 
pl.xlabel(r"$\omega_n$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Giw.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--------------------------------------VERTEX PLOTTING ------------------------------------------

print("Plotting vertex ...")

fdim = revert.shape[0]

if fdim <= shift:
    sys.exit("Error: Shift to large for vertex grid"); 

#---  Helper functions
def neg( w ):
    return fdim - w - 1

def check_bounds( w1, w2, w1p ):
    if ( w1 < 0 or w1 > fdim - 1 or w2 < 0 or w2 > fdim - 1 or w1p < 0 or w1p > fdim - 1 ):
        return False
    return True

def ReVert( w1, w2, w1p, i, j, k, l ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return revert[w1,w2,w1p,0,0,0,i,j,k,l]

def ImVert( w1, w2, w1p, i, j, k, l ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return imvert[w1,w2,w1p,0,0,0,i,j,k,l]

def ReVertUpDown( w1, w2, w1p ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return revert[w1, w2, w1p, 0, 0, 0, 0, 1, 0, 1]
    #return revert[w1, w2, w1p, 0, 0, 0, 0, 0, 0, 0]

def ImVertUpDown( w1, w2, w1p ):
    if ( not check_bounds( w1, w2, w1p ) ):
            return float('nan')
    return revert[w1, w2, w1p, 0, 0, 0, 0, 1, 0, 1]
    #return imvert[w1, w2, w1p, 0, 0, 0, 0, 0, 0, 0]

def ReVertUpUp( w1, w2, w1p ):
    return ReVert( w1, w2, w1p, 0, 1, 0, 1 ) - ReVert( w2, w1, w1p, 1, 0, 0, 1 )
    #return ReVert( w1, w2, w1p, 0, 0, 0, 0 ) - ReVert( w2, w1, w1p, 0, 0, 0, 0 )

def ImVertUpUp( w1, w2, w1p ):
    return ImVert( w1, w2, w1p, 0, 1, 0, 1 ) - ImVert( w2, w1, w1p, 1, 0, 0, 1 )
    #return ImVert( w1, w2, w1p, 0, 0, 0, 0 ) - ImVert( w2, w1, w1p, 0, 0, 0, 0 )


def plotVert( use_pl, zarr, string ):
    use_pl.set_aspect(1.0)
    pl.pcolormesh( vertgrid, vertgrid, np.ma.masked_where( np.isnan(zarr), zarr ) )
    pl.ylim([min(vertgrid),max(vertgrid)])
    pl.xlim([min(vertgrid),max(vertgrid)])
    use_pl.set_title( string , fontsize=10)
    pl.colorbar(shrink=0.6) 
    return

def plotUpUpVertRePP( use_pl ):
    title = RE + "\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{PP}-\omega_n,\omega_m)$"
    zarr = np.array([[ ReVertUpUp(n,shift+neg(n),m) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpDownVertRePP( use_pl ):
    title = RE + "\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{PP}-\omega_n,\omega_m)$"
    zarr = np.array([[ ReVertUpDown(n,shift+neg(n),m) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpUpVertRePH( use_pl ):
    title = RE + "\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n)$"
    zarr = np.array([[ ReVertUpUp(n,m+shift,n+shift) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpDownVertRePH( use_pl ):
    title = RE + "\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n)$"
    zarr = np.array([[ ReVertUpDown(n,m+shift,n+shift) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpUpVertReXPH( use_pl ):
    title = RE + "\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{XPH}+\omega_m,\omega_m)$"
    zarr = np.array([[ ReVertUpUp(n,m+shift,m) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpDownVertReXPH( use_pl ):
    title = RE + "\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{XPH}+\omega_m,\omega_m)$"
    zarr = np.array([[ ReVertUpDown(n,m+shift,m) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

#--- Plot

plotUpUpVertRePP( pl.subplot(2,3,1) )
pl.ylabel(r"$\omega_m$")
plotUpUpVertRePH( pl.subplot(2,3,2) )
plotUpUpVertReXPH( pl.subplot(2,3,3) )
plotUpDownVertRePP( pl.subplot(2,3,4) )
pl.xlabel(r"$\omega_n$")
pl.ylabel(r"$\omega_m$")
plotUpDownVertRePH( pl.subplot(2,3,5) )
pl.xlabel(r"$\omega_n$")
plotUpDownVertReXPH( pl.subplot(2,3,6) )
pl.xlabel(r"$\omega_n$")
pl.tight_layout()

#--- Save to file
pl.savefig("plots/Vert.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=2

#--- Plot Physical
plotUpUpVertRePP( pl.subplot(2,3,1) )
pl.ylabel(r"$\omega_m$")
plotUpUpVertRePH( pl.subplot(2,3,2) )
plotUpUpVertReXPH( pl.subplot(2,3,3) )

plotUpDownVertRePP( pl.subplot(2,3,4) )
pl.xlabel(r"$\omega_n$")
pl.ylabel(r"$\omega_m$")
plotUpDownVertRePH( pl.subplot(2,3,5) )
pl.xlabel(r"$\omega_n$")
plotUpDownVertReXPH( pl.subplot(2,3,6) )
pl.xlabel(r"$\omega_n$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Vert_shift.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

shift=0

