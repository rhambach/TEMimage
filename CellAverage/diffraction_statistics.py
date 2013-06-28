# tested influence of
#  - Nr:   no
#  - Nphi: no
#  - dp:   Bragg-Peak intensity increases with dp
#          seems to be due to some stray intensity
#          from Bragg Peak that is cut-off (reduced
#          if line-artefacts in FFT are reduced by
#          edge smoothin)
#  - binning: small influence (5%)
#
# Question: should we really use Integral of Bragg Peak?
#           or just the amplitude ?
#           How stable is the method ?
#           How does it relate to contrast measurements ?

import numpy as np
import matplotlib.pylab as plt
import wqplot.wq_stack as wq
import polartrafo 
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

def power_spectrum(img,scale=1,edge=0,verbosity=1):
  """
  Calculate Modulus of the Fourier Transform of an image.

    img   ... input image (array of shape NxM)
    scale ... (scalar or tuple) pixel size in x and y direction
    edge  ... (opt) width of edge which is smoothed before FFT
              to avoid artefacts due to jump at periodic border
    verbosity.(opt) quiet (0), debug (3)

  RETURNS absfft, qx, qy
    the origin is shifted to the center at O=(int(N/2), int(M/2))
  """
  # init parameters
  N,M = img.shape; 
  dx,dy = (scale,scale) if np.isscalar(scale) else scale;

  # mask edge
  if edge>0:
    assert edge < N/2 and edge < M/2  # edge should be small enough
    mx = np.zeros(N); my= np.zeros(M);
    mx[0:edge] = mx[-1:-edge-1:-1] = 1 - np.linspace(0,1,edge);
    my[0:edge] = my[-1:-edge-1:-1] = 1 - np.linspace(0,1,edge);
    mask=np.maximum(np.tile(mx,(M,1)).T,np.tile(my,(N,1)));
    I = np.mean(img);
    img = I*mask + img*(1-mask);  # blend between img and <img>

  # abs(fft)
  absfft = np.abs(np.fft.fftshift(np.fft.fft2(img)));
  qx     = np.fft.fftshift(np.fft.fftfreq(N))*2*np.pi/dx;
  qy     = np.fft.fftshift(np.fft.fftfreq(M))*2*np.pi/dy;

  # DEBUG test position of origin (N/2,M/2)
  assert abs(qx[1]-qx[0] - 2*np.pi/dx/N) < 1e-10
  assert qx[int(N/2)]== 0
  assert abs(qy[1]-qy[0] - 2*np.pi/dy/M) < 1e-10
  assert qy[int(M/2)]== 0
  assert absfft[int(N/2),int(M/2)] - np.mean(img)*N*M < 1e-10;
  if verbosity>3 and edge>0:
    plt.figure();
    plt.title("masked Image for FFT, edge width %d px"%edge)
    plt.imshow(img,aspect='auto',cmap=plt.cm.gray)

  return absfft, qx, qy

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
  filename="tests/TiO2_462eV.tif";
  img = tiff.imread(filename); # bin image to cut reciprocal radius
  #img = img[::2]+img[1::2];
  #img = img[:,::2]+img[:,1::2];# rebin by 2

  #img = img[1:]    # test for N<>M, N odd
  N,M = img.shape;  assert M==N;
  scale=0.015366*4            # nm/px
  scale*=10;                   #  A/px

  # perform FFT
  Ixy,qx,qy=power_spectrum(img,scale=scale,edge=0,verbosity=5);
  dqx=qx[1]-qx[0];
  dqy=qy[1]-qy[0];             # A-1/px

  # polar trafo
  rmin = 0.5;
  rmax = 5;                  # range for |q| [A-1]
  Irphi,r,phi = polartrafo.polar_distribution(Ixy,scale=(dqx,dqy),
                  Nphi=360,Nr=1000,rmin=rmin,rmax=rmax,verbosity=3);

  # plot polar distribution and fit peaks
  info={'xlabel': '|q|','ylabel':'polar angle',
        'xperchan': r[1]-r[0], 'yperchan': phi[1]-phi[0],
        'xoffset': r[0],'yoffset': phi[0],
        'xunits' : '1/A','yunits':'degree', 
        'desc':'polar distribution q*f(q,phi)'}
  fig1= wq.WQBrowser(Irphi,info, interpolation='none',aspect='auto');
  #x,y = find_peaks(Irphi,N=200);
  #fig1.axis.plot(r[x],phi[y], 'rx');
  
  dp     = 5
  bragg1 = 34 + np.asarray([-180,-90,0,90]); # define angles of bragg peaks
  bragg2 = bragg1+45

  def get_mask(bragg, c=(1,0,0,0.3)):
    mask   = np.zeros(len(phi),dtype=bool);
    for p0 in bragg:
      mask[np.abs(phi-p0)<dp] = True;
    # plot mask for reference
    color = np.zeros((len(phi),len(r),4));
    color[mask] = c; 
    fig1.axis.imshow(color, extent=fig1.imginfo['extent'],
                       interpolation='none',aspect='auto');
    return mask

  mask1  = get_mask(bragg1,c=(1,0,0,0.3));
  mask2  = get_mask(bragg2,c=(0,1,0,0.3));
  mask   = np.logical_or(mask1,mask2);

  # calculate average for background
  tot    = np.mean(Irphi,      axis=0);
  bg     = np.sum(Irphi[~mask],axis=0) / np.sum(~mask);  # average background spectrum
  sig1_bg= np.sum(Irphi[mask1],axis=0) / np.sum(mask1);  # average signal1+bg spectrum
  sig2_bg= np.sum(Irphi[mask2],axis=0) / np.sum(mask2);  # average signal2+bg spectrum

  plt.figure()
  plt.plot(r,tot,'k',label='average')
  plt.plot(r,bg, 'b',label='background');
  plt.plot(r,sig1_bg, 'r',label='first Bragg');
  plt.plot(r,sig2_bg, 'g',label='second Bragg');
  plt.ylim(0, np.max(sig1_bg[r>1]) + np.max(sig2_bg[r>1]));
  plt.legend();

  # statistics
  sig1   = ( sig1_bg - bg )  * np.sum(mask1);   # total intensity in 1. Bragg
  sig2   = ( sig2_bg - bg )  * np.sum(mask2);   # total intensity in 2. Bragg
  plt.figure();
  plt.plot(r,sig1)
  plt.plot(r,sig2)
  b1 = np.sum(sig1[np.logical_and(1.2<r, r<1.5)]);  
  b2 = np.sum(sig2[np.logical_and(1.8<r, r<2.1)]);  
  b3 = np.sum(sig1[np.logical_and(2.6<r, r<2.9)]);# problematic due to 4th Bragg in bg  

  print "# Statistics for Bragg Peak Intensity"
  print "#  image : '%s'" % filename
  print "#  dphi  : +/- %5.1f deg" % dp
  print "#  Bragg1  Bragg2  Bragg3  <Img>"
  print "   %f      %f      %f      %f   " % (b1,b2,b3,np.sum(Ixy))
  
  plt.show();
