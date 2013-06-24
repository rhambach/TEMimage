import numpy as np
import matplotlib.pylab as plt
import tifffile as tiff
import wqplot.wq_stack as wq
import img_filter as filt


def sum_unit_cell(img,a1,a2,N1,N2):
  y,x=np.ix_(np.fft.fftfreq(N)*N,np.fft.fftfreq(M)*M);
  shift1 = np.exp(-2j*np.pi/N*a1[0])**x * np.exp(-2j*np.pi/M*a1[1])**y
  shift2 = np.exp(-2j*np.pi/N*a2[0])**x * np.exp(-2j*np.pi/M*a2[1])**y

  shift=np.zeros((N,M),dtype=complex);
  for i in range(-N1/2,N1/2+1):  
    for j in range(-N2/2,N2/2+1):
      shift+= shift1**i * shift2**j
  shifted = np.fft.ifft2( np.fft.fft2(img) * shift );
  print "sum_unit_cell(): max of summed image, real=%5.3f, imag=%5.3f"\
                             %(np.max(shifted.real),np.max(shifted.imag));
  return shifted.real;


def plot_fft(img,b1,b2):
  # reciprocal lattice
  diff= np.fft.fftshift(np.fft.fft2(img));
  O   = np.asarray((N,M))/2.;   # center of image
  info={'xlabel': 'y','ylabel': 'x'}
  fig1= wq.WQBrowser(np.abs(diff),info,interpolation='none');


  b  = np.linalg.norm(b1);

  p = [(i*b1+j*b2+O) for i in range(-5,6) for j in range(-5,6)];
  fig1.axis.scatter(*zip(*p),marker='x',color='r');

  return fig1;

def dual_basis(b1,b2):
  # real-space lattice, a_i*b_j = N delta_ij
  a1 = np.asarray((-b2[1],b2[0]))/np.dot(b2,b2);  # perpendicular to b2
  a1*= N / np.dot(a1,b1);
  a2 = np.asarray((-b1[1],b1[0]))/np.dot(b1,b1); 
  a2*= N / np.dot(a2,b2);
  a  = np.linalg.norm(a1);
  print 'lattice constant  [A]:   ', a*scale*10;
  assert np.allclose((np.dot(a1,b1),np.dot(a2,b2)), N);
  assert np.allclose((np.dot(a2,b1),np.dot(a1,b2)), 0);
  return a1, a2;

def plot_unit_cells(img,a1,a2):
  info={'xlabel': 'y','ylabel': 'x'}
  info['xunits']=info['yunits']='nm';
  info['xperchan']=info['yperchan']=scale;
  fig2=wq.WQBrowser(img,info,interpolation='none');

  x = np.asarray((-N,N))*scale;
  y = x*a1[1]/a1[0];
  lines1 = [ z+n*a2[i]*scale for n in range(-50,50) for i,z in [(0,x),(1,y)]];
  y = x*a2[1]/a2[0];
  lines2 = [ z+n*a1[i]*scale for n in range(-50,50) for i,z in [(0,x),(1,y)]];

  fig2.axis.plot(*lines1,color='r');
  fig2.axis.plot(*lines2,color='r');
  return fig2;

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

# -- main ----------------------------------------
if __name__ == '__main__':
  import tifffile as tiff

  # load tiff image
  #filename="data/27_TiO2_0eV_ausschnitt_bin2.tif";
  filename="data/26_TiO2_463eV_ausschnitt_bin2.tif";
  imgstack = tiff.imread(filename);
  if len(imgstack.shape)==2: imgstack=imgstack[np.newaxis,:];
  img = imgstack[0];
  N,M = img.shape;
  d   = 11;        # number of unit-cells to sum
  scale=0.020237;  # nm/px
  assert M==N;

  # reciprocal unit vectors
  b1 = np.asarray([21.18415951,  32.51758871 ]);    # img7,8 # units [2pi/N]
  b2 = np.asarray([-b1[1],b1[0]]);     # units [2pi/M]
  fig1=plot_fft(img,b1,b2);

  # real space unit vectors
  a1,a2 = dual_basis(b1,b2);
  fig2=plot_unit_cells(img,a1,a2);
  
  # create averaged image by summing over unit cells
  for n,img in enumerate(imgstack):
    ucavg=sum_unit_cell(img,a1,a2,d,d);
    info={'xlabel': 'y','ylabel': 'x'}
    info['xunits']=info['yunits']='nm';
    info['xperchan']=info['yperchan']=scale;
    if n==0: fig3 =wq.WQBrowser(ucavg,info,interpolation='none');

    # output image
    outfile=filename.split('/')[-1].split('.tif')[0]+'_summed%02dx%02d_n%02d.tif'%(d,d,n);
    output_tif(outfile,ucavg,dtype='float32');  # unsigned float 32 bit (to read with DM3)
    plt.show();
