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
vert_mul=1

#--------------------------------------MANAGING FILES ------------------------------------------

def run(command):
    output = subprocess.check_output(command, shell=True)
    return output

fname = "full.h5"

f = h5py.File(fname, "r")

#--------------------------------------READ PARAMETERS FROM FILE ------------------------------------------

BETA=10.0
shift=0

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
imsig = f["/Sig/RE"]
resig = f["/Sig/IM"]
fdim = resig.shape[0]
#vertgrid = np.array(f["/Vert/fgrid"])
siggrid = np.arange(fdim)# np.array(f["/Giw/fgrid"])

#--- Plot function
def plotSig( i, j, use_pl, arr, title ):
    pl.plot( siggrid, arr[:,0,i,j], 'bx', ms=3, mew=0.2)
    #pl.xlim([1.2*min(vertgrid),1.2*max(vertgrid)])
    use_pl.set_title( title + "(" + str(i) + ","  + str(j) + ")$")
    return

#--- Plot
plotSig( 0, 0, pl.subplot(2,2,1), resig, RE + "\Sigma" ) 
plotSig( 0, 0, pl.subplot(2,2,2), imsig, IM + "\Sigma" ) 

plotSig( 1, 1, pl.subplot(2,2,3), resig, RE + "\Sigma" ) 
pl.xlabel(r"$\omega_n$")
plotSig( 1, 1, pl.subplot(2,2,4), imsig, IM + "\Sigma" ) 
pl.xlabel(r"$\omega_n$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Sig.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--------------------------------------G(iw) PLOTTING ------------------------------------------

print("Plotting Green function ...")


#--- Read
imGiw = f["/Giw/RE"]
reGiw = f["/Giw/IM"]
fdim = reGiw.shape[0]
Giwgrid = np.arange(fdim)# np.array(f["/Giw/fgrid"])

#--- Plot function
def plotGiw( i, j, use_pl, arr, title ):
    pl.plot( Giwgrid, arr[:,0,i,j], 'bx', ms=3, mew=0.2)
    #pl.xlim([1.2*min(vertgrid),1.2*max(vertgrid)])
    use_pl.set_title( title + str(i) + ","  + str(j) + ")$")
    return

#--- Plot
plotGiw( 0, 0, pl.subplot(2,2,1), reGiw, RE + "G(i\omega," ) 
plotGiw( 0, 0, pl.subplot(2,2,2), imGiw, IM + "G(i\omega," ) 

plotGiw( 1, 1, pl.subplot(2,2,3), reGiw, RE + "G(i\omega," ) 
pl.xlabel(r"$\omega_n$")
plotGiw( 1, 1, pl.subplot(2,2,4), imGiw, IM + "G(i\omega," ) 
pl.xlabel(r"$\omega_n$")

pl.tight_layout()


#--- Save to file
pl.savefig("plots/Giw.png", dpi=150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--------------------------------------VERTEX PLOTTING ------------------------------------------

print("Plotting vertex ...")


#--- Read
imvert = vert_mul * np.array(f["/Vert/RE"])
revert = vert_mul * np.array(f["/Vert/IM"])
fdim = revert.shape[0]
vertgrid = np.arange(fdim)# np.array(f["/Giw/fgrid"])

if fdim <= shift:
    sys.exit("Error: Shift to large for vertex grid"); 

#---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
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
    #return ReVert( w1, w2, w1p, 1, 0, 1, 0 )
    return ReVert( w1, w2, w1p, 0, 1, 0, 1 )

def ImVertUpDown( w1, w2, w1p ):
    #return ImVert( w1, w2, w1p, 1, 0, 1, 0 )
    return ImVert( w1, w2, w1p, 0, 1, 0, 1 )

def ReVertUpUp( w1, w2, w1p ):
    return ReVert( w1, w2, w1p, 0, 0, 0, 0 )

def ImVertUpUp( w1, w2, w1p ):
    return ImVert( w1, w2, w1p, 0, 0, 0, 0 )


def plotVert( use_pl, zarr, string ):
    use_pl.set_aspect(1.0)
    pl.pcolormesh( vertgrid, vertgrid, np.ma.masked_where( np.isnan(zarr), zarr ) )
    pl.ylim([min(vertgrid),max(vertgrid)])
    pl.xlim([min(vertgrid),max(vertgrid)])
    use_pl.set_title( string , fontsize=10)
    pl.colorbar(shrink=0.6) 
    return

def plotNambuVertRe( i, j, k, l, use_pl ):
    title = r"$\operatorname{Re}\gamma_2(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n," + str(i) + "," + str(j) + "," + str(k) + "," + str(l) + ")$" 
    zarr = np.array([[ ReVert(n,m,n,i,j,k,l) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title ) 
    return

def plotNambuVertIm( i, j, k, l, use_pl ):
    title = r"$\operatorname{Im}\gamma_2(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n," + str(i) + "," + str(j) + "," + str(k) + "," + str(l) + ")$" 
    zarr = np.array([[ ImVert(n,m,n,i,j,k,l) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpUpVertRePP( use_pl ):
    title = r"$\operatorname{Re}\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{PP}-\omega_n,\omega_m)$"
    zarr = np.array([[ ReVertUpUp(n,shift+neg(n),m) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpDownVertRePP( use_pl ):
    title = r"$\operatorname{Re}\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{PP}-\omega_n,\omega_m)$"
    zarr = np.array([[ ReVertUpDown(n,shift+neg(n),m) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpUpVertRePH( use_pl ):
    title = r"$\operatorname{Re}\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n)$"
    zarr = np.array([[ ReVertUpUp(n,m+shift,n+shift) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpDownVertRePH( use_pl ):
    title = r"$\operatorname{Re}\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n)$"
    zarr = np.array([[ ReVertUpDown(n,m+shift,n+shift) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpUpVertRePHX( use_pl ):
    title = r"$\operatorname{Re}\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{PHX}+\omega_m,\omega_m)$"
    zarr = np.array([[ ReVertUpUp(n,m+shift,m) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return

def plotUpDownVertRePHX( use_pl ):
    title = r"$\operatorname{Re}\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{PHX}+\omega_m,\omega_m)$"
    zarr = np.array([[ ReVertUpDown(n,m+shift,m) for n in range(fdim)] for m in range(fdim)])
    plotVert( use_pl, zarr, title )
    return


#pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega_{\rm PP}=\Omega_{\rm PH}=\Omega_{\rm xPH}=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\gamma_{2}(\omega_1,\omega_2,\omega_1')$")


#--- Plot Physical
plotUpUpVertRePP( pl.subplot(2,3,1) )
pl.ylabel(r"$\omega_m$")
plotUpUpVertRePH( pl.subplot(2,3,2) )
plotUpUpVertRePHX( pl.subplot(2,3,3) )

plotUpDownVertRePP( pl.subplot(2,3,4) )
pl.xlabel(r"$\omega_n$")
pl.ylabel(r"$\omega_m$")
plotUpDownVertRePH( pl.subplot(2,3,5) )
pl.xlabel(r"$\omega_n$")
plotUpDownVertRePHX( pl.subplot(2,3,6) )
pl.xlabel(r"$\omega_n$")

pl.tight_layout()

#--- Save to file
pl.savefig("plots/Vert.png", dpi = 150)
pl.figure(dpi=100) # Reset dpi to default
pl.clf()

#--------------------------------------G2 PLOTTING ------------------------------------------

#print("Plotting G2 ...")


##--- Read
#revert = vert_mul * np.array(f["/G2/RE"])
#imvert = vert_mul * np.array(f["/G2/IM"])
#fdim = revert.shape[0]
#vertgrid = np.arange(fdim)# np.array(f["/Giw/fgrid"])

#if fdim <= shift:
    #sys.exit("Error: Shift to large for vertex grid"); 

##---  Helper functions (include translation from  nambu to physical CAUTION : USING TIME REVERSAL SYMM) 
#def neg( w ):
    #return fdim - w - 1

#def check_bounds( w1, w2, w1p ):
    #if ( w1 < 0 or w1 > fdim - 1 or w2 < 0 or w2 > fdim - 1 or w1p < 0 or w1p > fdim - 1 ):
        #return False
    #return True

#def ReVert( w1, w2, w1p, i, j, k, l ):
    #if ( not check_bounds( w1, w2, w1p ) ):
            #return float('nan')
    #return revert[w1,w2,w1p,0,0,0,i,j,k,l]

#def ImVert( w1, w2, w1p, i, j, k, l ):
    #if ( not check_bounds( w1, w2, w1p ) ):
            #return float('nan')
    #return imvert[w1,w2,w1p,0,0,0,i,j,k,l]

#def ReVertUpDown( w1, w2, w1p ):
    #return ReVert( w1, w2, w1p, 1, 0, 1, 0 )
    ##return ReVert( w1, w2, w1p, 0, 1, 0, 1 )

#def ImVertUpDown( w1, w2, w1p ):
    #return ImVert( w1, w2, w1p, 1, 0, 1, 0 )
    ##return ImVert( w1, w2, w1p, 0, 1, 0, 1 )

#def ReVertUpUp( w1, w2, w1p ):
    #return ReVert( w1, w2, w1p, 0, 0, 0, 0 )

#def ImVertUpUp( w1, w2, w1p ):
    #return ImVert( w1, w2, w1p, 0, 0, 0, 0 )


#def plotVert( use_pl, zarr, string ):
    #use_pl.set_aspect(1.0)
    #pl.pcolormesh( vertgrid, vertgrid, np.ma.masked_where( np.isnan(zarr), zarr ) )
    #pl.ylim([min(vertgrid),max(vertgrid)])
    #pl.xlim([min(vertgrid),max(vertgrid)])
    #use_pl.set_title( string , fontsize=10)
    #pl.colorbar(shrink=0.6) 
    #return

#def plotNambuVertRe( i, j, k, l, use_pl ):
    #title = r"$\operatorname{Re}\gamma_2(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n," + str(i) + "," + str(j) + "," + str(k) + "," + str(l) + ")$" 
    #zarr = np.array([[ ReVert(n,m,n,i,j,k,l) for n in range(fdim)] for m in range(fdim)])
    #plotVert( use_pl, zarr, title ) 
    #return

#def plotNambuVertIm( i, j, k, l, use_pl ):
    #title = r"$\operatorname{Im}\gamma_2(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n," + str(i) + "," + str(j) + "," + str(k) + "," + str(l) + ")$" 
    #zarr = np.array([[ ImVert(n,m,n,i,j,k,l) for n in range(fdim)] for m in range(fdim)])
    #plotVert( use_pl, zarr, title )
    #return

#def plotUpUpVertRePP( use_pl ):
    #title = r"$\operatorname{Re}\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{PP}-\omega_n,\omega_m)$"
    #zarr = np.array([[ ReVertUpUp(n,shift+neg(n),m) for n in range(fdim)] for m in range(fdim)])
    #plotVert( use_pl, zarr, title )
    #return

#def plotUpDownVertRePP( use_pl ):
    #title = r"$\operatorname{Re}\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{PP}-\omega_n,\omega_m)$"
    #zarr = np.array([[ ReVertUpDown(n,shift+neg(n),m) for n in range(fdim)] for m in range(fdim)])
    #plotVert( use_pl, zarr, title )
    #return

#def plotUpUpVertRePH( use_pl ):
    #title = r"$\operatorname{Re}\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n)$"
    #zarr = np.array([[ ReVertUpUp(n,m+shift,n+shift) for n in range(fdim)] for m in range(fdim)])
    #plotVert( use_pl, zarr, title )
    #return

#def plotUpDownVertRePH( use_pl ):
    #title = r"$\operatorname{Re}\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{PH}+\omega_m,\Omega_{PH}+\omega_n)$"
    #zarr = np.array([[ ReVertUpDown(n,m+shift,n+shift) for n in range(fdim)] for m in range(fdim)])
    #plotVert( use_pl, zarr, title )
    #return

#def plotUpUpVertRePHX( use_pl ):
    #title = r"$\operatorname{Re}\gamma_{2,\uparrow\uparrow}(\omega_n,\Omega_{PHX}+\omega_m,\omega_m)$"
    #zarr = np.array([[ ReVertUpUp(n,m+shift,m) for n in range(fdim)] for m in range(fdim)])
    #plotVert( use_pl, zarr, title )
    #return

#def plotUpDownVertRePHX( use_pl ):
    #title = r"$\operatorname{Re}\gamma_{2,\uparrow\downarrow}(\omega_n,\Omega_{PHX}+\omega_m,\omega_m)$"
    #zarr = np.array([[ ReVertUpDown(n,m+shift,m) for n in range(fdim)] for m in range(fdim)])
    #plotVert( use_pl, zarr, title )
    #return


##pl.suptitle(r"$U=$" + str(UINT) + r"     $\beta=$" + str(BETA) + r"     $\epsilon=$" + str(EPS) +  r"     $\Omega_{\rm PP}=\Omega_{\rm PH}=\Omega_{\rm xPH}=$" + str(shift) + r"$*2\pi/\beta$" + r"     Notation: $\gamma_{2}(\omega_1,\omega_2,\omega_1')$")


##--- Plot Physical
#plotUpUpVertRePP( pl.subplot(2,3,1) )
#pl.ylabel(r"$\omega_m$")
#plotUpUpVertRePH( pl.subplot(2,3,2) )
#plotUpUpVertRePHX( pl.subplot(2,3,3) )

#plotUpDownVertRePP( pl.subplot(2,3,4) )
#pl.xlabel(r"$\omega_n$")
#pl.ylabel(r"$\omega_m$")
#plotUpDownVertRePH( pl.subplot(2,3,5) )
#pl.xlabel(r"$\omega_n$")
#plotUpDownVertRePHX( pl.subplot(2,3,6) )
#pl.xlabel(r"$\omega_n$")

#pl.tight_layout()

##--- Save to file
#pl.savefig("plots/G2.png", dpi = 150)
#pl.figure(dpi=100) # Reset dpi to default
#pl.clf()

