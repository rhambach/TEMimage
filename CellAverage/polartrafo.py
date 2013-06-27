 ######################################
# Tests for pyFAI
#  - works on discretized probability densities, i.e., input is f(x,y)
#    output is f(r,phi) or f(r) for integrate2d/1d
#  - connection with image intensity (binned probability density):
#    cartesian coords
#     Ixy =  Dx Dy f(x,y)    (Image value at point x,y, pixel size Dx,Dy)
#     I   =  Int dx dy f(x,y) = Sum Ixy  (total image intensity)
#    polar coords
#     Irphi= Dr Dphi r f(r,phi)  
#     I   =  Int dr dphi r f(r,phi) = Sum Irphi
#  - different units are used in pixel size for __init__ (m) and integrate?d (mm)
#  - units for angle seem to be not always consistently degree or rad internally
#  - integrate1d and integrate2d are equivalent as long as we do not 
#    exceed the image size r<L/sqrt(2)
########################################

import pyFAI;
import numpy as np
import matplotlib.pylab as plt
import wqplot.wq_stack as wq

def polar_distribution(Ixy, O=None, scale=1, rmin=0, rmax=None, Nr=100, Nphi=360, verbosity=1):
  """ 
  wrapper for transformation from cartesian to polar coordinates 
  uses pyFAI.integrate2d()

   Ixy   ...  image in cartesian coordinates
   O     ...  (tuple) origin of coordinate system [px]
                      (default: center of image)
   scale ...  (scalar or tuple) pixel size in x and y direction
   rmin
   rmax  ...  range for radius in output array
   Nr    ...  number of pixels along radial direction
   Nphi  ...  number of pixels along angular direction
   verbosity. 0: silent, 1: normal, 3: debug

  RETURNS
   Irphi, r, phi ... arrays containing transformed image and 
                     the corresponding axes
  """ 
  # init parameters
  N,M = Ixy.shape; 
  dx,dy = (scale,scale) if np.isscalar(scale) else scale;
  if O is None:       # default O: center of img 
    O = np.asarray((int(N/2),int(M/2)));
  if rmax is None:    # default rmax: largest circle in img
    rmax = min(N/2.*dx,M/2.*dy);

  # setup polar trafo
  # NOTE: pixel size in __init__ is given in meter, while
  #       all other units are returned in mm
  dr   = (rmax-rmin)/float(Nr);
  dphi = 2*np.pi/float(Nphi);
  ai = pyFAI.AzimuthalIntegrator(pixel1=dx/1000.,pixel2=dy/1000.,
          poni1=O[0]*dx/1000.,poni2=O[1]*dy/1000.); # pixel size in m

  # wrapper to calculate polar distribution f(r,phi)
  fxy   = Ixy/dx/dy;           # density distribution f(x,y)
  frphi,r,phi=ai.integrate2d(fxy,Nr,Nphi,radial_range=(rmin,rmax),
          unit=pyFAI.units.R); # units in mm
  Irphi = frphi*r*dphi*dr      # binned intensity distribution I[r,phi]
 
  # DEBUG
  if verbosity > 2:
    print "Test total intensity of image Ixy: ", np.sum(Ixy)
    print "      and transformed image Irphi: ", np.sum(Irphi)
  assert np.allclose(r[1:]-r[:-1],dr);       # test step sizes
  assert np.allclose(phi[1:]-phi[:-1],360/float(Nphi),rtol=0.01);

  return Irphi, r, phi



def radial_distribution(Ixy, O=None, scale=1, rmin=0, rmax=None, Nr=100, verbosity=1):
  """ 
  wrapper to determine radial distribution f(r)*dr for an image
  uses pyFAI.integrate1d(). The total intensity is given by
   
   I = Sum(Ixy) = Int dr r f(r) = Sum r frdr 

   Ixy   ...  image in cartesian coordinates
   O     ...  (tuple) origin of coordinate system [px]
                      (default: center of image)
   scale ...  (scalar or tuple) pixel size in x and y direction
   rmin
   rmax  ...  range for radius in output array
   Nr    ...  number of pixels along radial direction
   verbosity. 0: silent, 1: normal, 3: debug

  RETURNS
   fr, r ... arrays containing the radial distribution transformed image and 
               the corresponding axes
  """ 
  # init parameters
  N,M = Ixy.shape; 
  dx,dy = (scale,scale) if np.isscalar(scale) else scale;
  if O is None:       # default O: center of img 
    O = np.asarray((int(N/2),int(M/2)));
  if rmax is None:    # default rmax: largest circle in img
    rmax = min(N/2.*dx,M/2.*dy);

  # setup polar trafo
  # NOTE: pixel size in __init__ is given in meter, while
  #       all other units are returned in mm
  dr   = (rmax-rmin)/float(Nr);
  dphi = 2*np.pi/float(Nphi);
  ai = pyFAI.AzimuthalIntegrator(pixel1=dx/1000.,pixel2=dy/1000.,
          poni1=O[0]*dx/1000.,poni2=O[1]*dy/1000.); # pixel size in m

  # radial distribution f(r) = dr * Int dphi f(r,phi)
  fxy   = Ixy/dx/dy;           # density distribution f(x,y)
  r,fr  = ai.integrate1d(fxy,Nr,radial_range=(rmin,rmax),
          unit=pyFAI.units.R); # in mm
  frdr  = fr* 2*np.pi*dr;      # correct normalisation
 
  # DEBUG
  if verbosity > 2:
    print "Test total intensity of image Ixy: ", np.sum(Ixy)
    print " and the radial distribution frdr: ", np.sum(frdr*r)
  assert np.allclose(r[1:]-r[:-1],dr);       # test step sizes
 
  return frdr, r


# -- main ----------------------------------------
if __name__ == '__main__':
  import tifffile as tiff

  # load tiff image
  filename="tests/TiO2_0eV.tif";
  img = tiff.imread(filename); # bin image to cut reciprocal radius
  N,M = img.shape;  assert M==N;
  Scale=0.015366*4;            # nm/px
  Dx=Dy= 2*np.pi/Scale/N/10;   # bin width for cartesian coord [A-1/px]
  edge = 20                    # width of edge to be smoothed

  # reciprocal lattice
  Ixy = np.abs(np.fft.fftshift(np.fft.fft2(img)));
  info={'xlabel': 'y','ylabel': 'x',
        'desc'  : 'Image in cartesian coordinates'}
  fig1= wq.WQBrowser(Ixy,info,interpolation='none');
  fig1.set_style("log");

  # polar trafo
  rmin = 0;
  rmax = None;                  # range for |q| [A-1]
  Nr   = 100;                # number of bins along radius
  Nphi = 360;                # number of bins along angle
  Irphi, r2, phi = polar_distribution(Ixy,scale=(Dx,Dy),Nr=Nr,Nphi=Nphi,rmin=rmin,rmax=rmax);

  # radial distribution
  frdr, r1 = radial_distribution(Ixy,scale=Dx,Nr=Nr,rmin=rmin,rmax=rmax);

  # comparision polar and radial trafo
  plt.figure();             
  plt.title('radial distribution');
  plt.xlabel('|q| [1/A]'); plt.ylabel('q*f(q)')
  plt.plot(r1,frdr*r1,label="f(r) r dr");
  plt.plot(r2,np.sum(Irphi,axis=0),label="Int r f(r,phi) dphi dr");
  plt.legend();
  print "Test total intensity of image: ", np.sum(Ixy)
  print "            - sum integrate2d: ", np.sum(Irphi)
  print "            - sum integrate1d: ", np.sum(frdr*r1)

  # plot polar distribution and fit peaks
  info={'xlabel': '|q|','ylabel':'polar angle',
        'xperchan': r2[1]-r2[0], 'yperchan': phi[1]-phi[0],
        'xoffset': r2[0],'yoffset': phi[0],
        'xunits' : '1/A','yunits':'degree', 
        'desc':'polar distribution q*f(q,phi)'}
  fig2= wq.WQBrowser(Irphi,info, interpolation='none',aspect='auto');
  plt.show();
