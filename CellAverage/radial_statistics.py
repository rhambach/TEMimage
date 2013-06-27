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

import numpy as np
import matplotlib.pylab as plt
import tifffile as tiff
import wqplot.wq_stack as wq
import img_filter as filt
import pyFAI;
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

def gauss2D(x,y, x0,y0, A, fwhm, bg):
  """
   normalised 2D gaussian function 
   A    ... integral
   x0,y0... central position
   fwhm ... full widht half maximum of marginal probabilities
   bg   ... background
   [http://mathworld.wolfram.com/GaussianFunction.html]
  """
  sigma = fwhm / 2.354820045;  # fwhm / 2*sqrt(2 ln(2))  
  return bg + A/(2*np.pi*sigma**2) * \
          np.exp(-((x-x0)**2+(y-y0)**2)/(2*sigma**2));

def output_tif(outfile,img,int8=False,dtype=None):
  import tifffile as tiff
  if int8:
    I = img.copy();
    I = (I - I.min())/(I.max()-I.min())*255;  # normalize to [0,256);
    tiff.imsave(outfile,I.astype('int8'));
  elif dtype is not None:
    I = img.astype(dtype);
    print "output_tif(): limits in img: ", I.min(), I.max();
    tiff.imsave(outfile,I);
  else:
    tiff.imsave(outfile,img);
plt.show(); 

def find_objects(data, size=5, thresh=None, N=20, verbosity=1):

  # maximum step height in neighborhood of size 'size'
  data_max = filters.maximum_filter(data, size)
  data_min = filters.minimum_filter(data, size)
  diff = data_max - data_min;

  # determine threshold which gives N pixels
  if thresh is None:  thresh=np.sort(diff.flat)[-N];

  # create mask for image (maximum with high step size)
  maxima = (data == data_max)
  maxima[diff <= thresh] = False
  
  # find connected objects (pixel agglomerates)
  labeled, num_objects = ndimage.label(maxima)
  slices = ndimage.find_objects(labeled)

  # DEBUG
  if verbosity>2:
    plt.figure()
    plt.imshow(np.log(1+np.abs(data)),interpolation='nearest');
    # create overlay with red color (r,g,b,alpha) at each image point 
    overlay = np.zeros(data.shape+(4,));
    overlay[maxima] = (1,1,0,1); # set opacity according to maxima
    plt.imshow(overlay,interpolation='nearest');

  return slices;

def find_peaks(data, size=5, fitradius=None, thresh=None, N=20,verbosity=1):
  # iterate to find really the number N ?
  slices = find_objects(data,size,thresh=thresh,N=N,verbosity=verbosity-10);  

  # iterate over peaks
  pos = [];
  for dy,dx in slices:

    # initial guess for peak parameters
    x = (dx.start + dx.stop - 1)/2
    y = (dy.start + dy.stop - 1)/2    
    
    # refined fit
    if fitradius is not None:
      X,Y=np.ogrid[0:2*fitradius+1,0:2*fitradius+1];
      X-=fitradius+x; Y-=fitradius+y;
      region    = data[X,Y],# (periodic boundary conditions)
      print x,y, fitradius, region.shape;
      plt.figure();
      plt.imshow(region);
      plt.show();
      param0    = (0,0,np.sum(data),1,0);
      residuals = lambda param,data:  gauss(X,Y,*param);

    pos.append((x,y));

  return np.transpose(pos);

def refine_fit(data,peaks,model,p0):
  for x,y in peaks:
    # fit model function
    pass

# -- main ----------------------------------------
if __name__ == '__main__':
  import tifffile as tiff

  # load tiff image
  filename="tests/TiO2_0eV.tif";
  img = tiff.imread(filename); # bin image to cut reciprocal radius
  N,M = img.shape;  assert M==N;
  scale=0.015366*4;            # nm/px
  dx=dy= 2*np.pi/scale/N/10;   # bin width for cartesian coord [A-1/px]

  # reciprocal lattice
  Ixy= np.abs(np.fft.fftshift(np.fft.fft2(img)));
  fxy= Ixy/dx/dy;              # density distribution f(x,y)
  I  = np.sum(Ixy)             # total intensity
  O  = np.asarray((N,M))/2.;   # center of image [px]
  info={'xlabel': 'y','ylabel': 'x'}
  fig1= wq.WQBrowser(Ixy,info,interpolation='none');
  print O

  # polar trafo
  rmin = 0;
  rmax = 5;                 # range for |q| [A-1]
  Nr   = 100;                # number of bins along radius
  Nphi = 360;                # number of bins along angle
  dr   = (rmax-rmin)/float(Nr);
  dphi = 2*np.pi/float(Nphi);
  
  # create PolarTrafo object and set pixel size
  # NOTE: pixel size in __init__ is given in meter, while
  #       all other units are returned in mm
  pw=dx/1000;
  ai = pyFAI.AzimuthalIntegrator(pixel1=pw,pixel2=pw,poni1=O[0]*pw,poni2=O[1]*pw); # pixel size in m

  # radial distribution f(r)/(2*pi) = Int dphi f(r,phi) / (2*pi)
  r,fr = ai.integrate1d(fxy,Nr,radial_range=(rmin,rmax),unit=pyFAI.units.R); # in mm
  plt.figure();
  plt.title('radial distribution');
  plt.xlabel('|q| [1/A]'); plt.ylabel('q*f(q)')
  plt.plot(r,fr*r*2*np.pi);

  # polar distribution f(r,phi)
  frphi,r,phi=ai.integrate2d(fxy,Nr,Nphi,radial_range=(rmin,rmax),unit=pyFAI.units.R); # mm
  Irphi = frphi*r*dphi*dr
  plt.plot(r,np.sum(frphi,axis=0)*dphi*r);   # comparison with radial distribution
   
  print "Test total intensity of image: ", I # test total intensity
  print "            - sum integrate2d: ", np.sum(Irphi)
  print "            - sum integrate1d: ", np.sum(fr*r)*dr*2*np.pi

  assert np.allclose(r[1:]-r[:-1],dr);       # test step sizes
  assert np.allclose(phi[1:]-phi[:-1],360/float(Nphi),rtol=0.01);

  # plot polar distribution and fit peaks
  info={'xlabel': '|q|','ylabel':'polar angle',
        'xperchan': r[1]-r[0], 'yperchan': phi[1]-phi[0],
        'xoffset': r[0],'yoffset': phi[0],
        'xunits' : '1/A','yunits':'degree', 
        'desc':'polar distribution q*f(q,phi)'}
  x,y = find_peaks(Irphi,N=2000);

  fig2= wq.WQBrowser(Irphi,info, interpolation='none',aspect='auto');
  fig2.axis.plot(r[x],phi[y], 'rx')

  plt.show();
