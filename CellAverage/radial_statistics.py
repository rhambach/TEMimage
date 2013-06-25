import numpy as np
import matplotlib.pylab as plt
import tifffile as tiff
import wqplot.wq_stack as wq
import img_filter as filt
import pyFAI;
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

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

def find_peaks(data, neighborhood_size=5, threshold=5e4):

  data_max = filters.maximum_filter(data, neighborhood_size)
  maxima = (data == data_max)
  data_min = filters.minimum_filter(data, neighborhood_size)
  diff = ((data_max - data_min) > threshold)
  maxima[diff == 0] = 0
  plt.figure()
  plt.imshow(-maxima,interpolation='nearest')
 
  labeled, num_objects = ndimage.label(maxima)
  slices = ndimage.find_objects(labeled)
  print len(slices);
  x, y = [], []
  for dy,dx in slices:
    x_center = (dx.start + dx.stop - 1)/2
    x.append(x_center)
    y_center = (dy.start + dy.stop - 1)/2    
    y.append(y_center)
  return x,y

def refine_fit(data,peaks,model,p0):
  for x,y in peaks:
    # fit model function
    pass

# -- main ----------------------------------------
if __name__ == '__main__':
  import tifffile as tiff

  # load tiff image
  filename="tests/TiO2_0eV.tif";
  img = tiff.imread(filename);
  N,M = img.shape;
  scale=0.015366*4;  # nm/px
  assert M==N;

  # reciprocal lattice
  diff= np.abs(np.fft.fftshift(np.fft.fft2(img)));
  O   = np.asarray((N,M))/2.;   # center of image
  info={'xlabel': 'y','ylabel': 'x'}
  fig1= wq.WQBrowser(diff,info,interpolation='none');
  print O

  # polar trafo
  dq = 2*np.pi/scale/N/10; # in 1/A/px
  dq/= 1000;            # take into account difference in units in pyFAI
                        # pixel size in m, results in mm
  ai = pyFAI.AzimuthalIntegrator(pixel1=dq,pixel2=dq,poni1=O[0]*dq,poni2=O[1]*dq); # pixel size in m
  out= ai.integrate1d(diff,1000,radial_range=(0,5),unit=pyFAI.units.R); # in mm
  plt.figure();
  plt.plot(*out);

  rphi,r,phi=ai.integrate2d(diff,300,radial_range=(0,5),unit=pyFAI.units.R); # mm
  dr=r[1]-r[0];
  dphi=phi[1]-phi[0];
  assert np.allclose(r[1:]-r[:-1],dr);
  assert np.allclose(phi[1:]-phi[:-1],dphi);
  info={'xlabel': 'x','ylabel':'y',
        'xperchan': dr, 'yperchan': dphi,
        'xoffset': r[0],'yoffset': phi[0],
        'xunits' : '1/A','yunits':'degree'}
  x,y = find_peaks(rphi);
  fig2= wq.WQBrowser(rphi,info,interpolation='none',aspect='auto');
  fig2.axis.plot(r[x],phi[y], 'rx')

  plt.show();
