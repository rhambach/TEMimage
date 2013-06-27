import numpy as np
import matplotlib.pylab as plt
import wqplot.wq_stack as wq
import polartrafo 
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

def power_spectrum(img,scale=1,edge=0):
  """
  Calculate Modulus of the Fourier Transform of an image.

    img   ... input image (array of shape NxM)
    scale ... (scalar or tuple) pixel size in x and y direction
    edge  ... (opt) width of edge which is smoothed before FFT
              to avoid artefacts due to jump at periodic border

  RETURNS absfft, qx, qy
    the origin is shifted to the center at O=(int(N/2), int(M/2))
  """
  # init parameters
  N,M = img.shape; 
  dx,dy = (scale,scale) if np.isscalar(scale) else scale;

  # mask edge
  if edge>0:
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

  assert abs(qx[1]-qx[0] - 2*np.pi/dx/N) < 1e-10
  assert qx[int(N/2)]== 0
  assert abs(qy[1]-qy[0] - 2*np.pi/dy/M) < 1e-10
  assert qy[int(M/2)]== 0

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
  img = tiff.imread(filename); # bin image to cut reciprocal radius
  N,M = img.shape;  assert M==N;
  scale=0.015366*4;            # nm/px
  scale*=10;                   #  A/px

  # perform FFT
  Ixy,qx,qy=power_spectrum(img,scale=scale,edge=20);
  dqx=qx[1]-qx[0];
  dqy=qy[1]-qy[0];             # A-1/px

  # polar trafo
  rmin = 0;
  rmax = 5;                 # range for |q| [A-1]
  Irphi,r,phi = polartrafo.polar_distribution(Ixy,scale=(dqx,dqy),
                  Nphi=360,Nr=100,rmin=rmin,rmax=rmax,verbosity=3);

  # plot polar distribution and fit peaks
  info={'xlabel': '|q|','ylabel':'polar angle',
        'xperchan': r[1]-r[0], 'yperchan': phi[1]-phi[0],
        'xoffset': r[0],'yoffset': phi[0],
        'xunits' : '1/A','yunits':'degree', 
        'desc':'polar distribution q*f(q,phi)'}
  fig1= wq.WQBrowser(Irphi,info, interpolation='none',aspect='auto');
  x,y = find_peaks(Irphi,N=2000);
  fig1.axis.plot(r[x],phi[y], 'rx')

  plt.show();
