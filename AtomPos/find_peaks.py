"""
  Interactive fitting of peaks in noisy 2D images.

  Copyright (c) 2013, rhambach. 
    This file is part of the FitAtomPos package and released
    under the MIT-Licence. See LICENCE file for details.
"""
import numpy as np
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt
import point_browser as pb
from   matplotlib.widgets import Button, RadioButtons, Slider

class FindCenters(pb.PointBrowser):
  """
  Semi-automatic fitting of bright points in TEM image
    dragging ... (opt) if True, dragging is allowed for sliders
  """
  def __init__(self, image, dragging=True, **kwargs):

    # init PointBrowser
    super(FindCenters,self).__init__(image,[[None,None]],**kwargs);
    self.axis.set_title('FitHexagonCenters: %s' % self.imginfo['desc']);
    self.fig.subplots_adjust(bottom=0.2);  # space for sliders

    # add slider for neighborhood size
    self.nbhd_size = 5;        # neighborhood size
    axNbhd = self.fig.add_axes([0.2, 0.05, 0.1, 0.04]);
    self.sNbhd = Slider(axNbhd,'neighbors ',2,20,valinit=self.nbhd_size,\
                                valfmt=' (%d)',dragging=dragging);
    self.sNbhd.on_changed(self.ChangeNeighborhood);

    # add slider for number of points
    self.num_points = 1e6;     # number of local maxima to find
    axNum  = self.fig.add_axes([0.45, 0.05, 0.3, 0.04]);
    self.sNum = Slider(axNum, 'points ',0,100,valinit=self.num_points,\
                               valfmt=' (%d)',dragging=dragging);
    self.sNum.on_changed(self.ChangeMaxPoints);

    # add buttons
    axRefine = self.fig.add_axes([0.85, 0.7, 0.1, 0.04]);
    self.bRefine = Button(axRefine,'Refine');
    self.bRefine.on_clicked(self.RefineCenters);

    # initial calculation of local maximas
    self.ChangeNeighborhood(self.nbhd_size);


  def ChangeNeighborhood(self,val):
    #print "ChangeNeighborhood"
    self.nbhd_size = int(val);
    # run initial peak fit
    maxima,diff = self.find_local_maxima(self.image,self.nbhd_size);
    # update max-number of points in points slider
    self.sNum.valmax = Nmax = np.sum(maxima);           # number of local max
    self.sNum.ax.set_xlim((self.sNum.valmin, Nmax));    # rescale slider
    self.sNum.set_val(min(self.num_points,Nmax));       # update value (calls ChangeMaxPoints)
 

  def ChangeMaxPoints(self,val):
    #print "ChangeMaxPoints()";
    self.num_points = int(val);
    self.points = self.refine_local_maxima(self.num_points);
    self._update_points();
    

  def RefineCenters(self,event):
    " refine positions by fitting 2D Gaussian in neighborhood of local max "
    from  scipy.optimize import leastsq
    from  sys import stdout

    #print "Refine()";
    NN    = self.nbhd_size;
    Nx,Ny = self.image.shape;
    dx,dy = np.mgrid[-NN:NN+1,-NN:NN+1];

    # refine each point separately
    self.points = self.points.astype(float); # allow subpixel precision
    for ip in range(len(self.points)):
      P   = self.points[ip];
      x,y = np.round(P);

      # get neighborhood (skip border)
      xmin,xmax = dx[[0,-1],0]+x;  # first and last element in dx
      ymin,ymax = dy[0,[0,-1]]+y;  #       "                   dy
      if xmin<0 or ymin<0 or xmax>=Nx or ymax>=Ny: continue
      nbhd = self.image[xmin:xmax+1,ymin:ymax+1];
      assert nbhd.shape == (2*NN+1,2*NN+1)

      # calculate center of mass
      def gauss(x0,y0,A,B,fwhm):
        return A*np.exp( - ((dx-x0)**2+(dy-y0)**2) / fwhm**2) + B;
      p0     = (0.,0.,self.image[tuple(P)],0.,NN/2);          # initial guess
      residuals = lambda param: (nbhd - gauss(*param)).flat;  # residuals
      p,ierr = leastsq(lambda p: (nbhd - gauss(*p)).flat, p0);# least-squares fit
      self.points[ip] = (x+p[0],y+p[1]);             # correct position of point

      # DEBUG: plot fits for each point
      if self.verbosity > 0:
        print "Refining Points...  %d %%\r" % (100*ip/len(self.points-1)),
      if self.verbosity > 3:
        print "IN:  ",p0
        print "OUT: ",p
      if self.verbosity > 10:
        plt.figure();
        ix = nbhd.shape[0]/2;
        plt.plot(dy[ix],nbhd[ix],  'k',label='image');
        plt.plot(dy[ix],gauss(*p0)[ix],'g',label='first guess');
        plt.plot(dy[ix],gauss(*p)[ix], 'r',label='final fit');
        plt.plot(dx[:,ix],nbhd[:,ix],      'k--');
        plt.plot(dx[:,ix],gauss(*p0)[:,ix],'g--');
        plt.plot(dx[:,ix],gauss(*p)[:,ix], 'r--');
        plt.legend();
        plt.show();
    if self.verbosity > 0:  print "Refining Points. Finished.";
    stdout.flush();

    self._update_points();


  def find_local_maxima(self, data, neighborhood_size):
    """ 
     find local maxima within neighborhood 
      idea from http://stackoverflow.com/questions/9111711
      (get-coordinates-of-local-maxima-in-2d-array-above-certain-value)
    """

    # find local maxima in image (width specified by neighborhood_size)
    data_max = filters.maximum_filter(data,neighborhood_size);
    maxima   = (data == data_max);
    assert np.sum(maxima) > 0;        # we should always find local maxima
  
    # remove connected pixels (plateaus)
    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
    maxima *= 0;
    for dx,dy in slices:
      maxima[(dx.start+dx.stop-1)/2, (dy.start+dy.stop-1)/2] = 1

    # calculate difference between local maxima and lowest 
    # pixel in neighborhood (will be used in select_local_maxima)
    data_min = filters.minimum_filter(data,neighborhood_size);
    diff     = data_max - data_min;
    self._maxima = maxima;
    self._diff   = diff;

    return maxima,diff

  def refine_local_maxima(self,N):
    " select highest N local maxima using thresholding "

    maxima = self._maxima;  diff = self._diff;

    # select highest local maxima using thresholding
    if np.sum(maxima) > N:
      # calc treshold from sorted list of differences for local maxima
      thresh = np.sort(diff[maxima].flat)[-N];
      # keep only maxima with diff>thresh
      maxima = np.logical_and(maxima, diff>thresh);  

    # TODO: refine fit by local 2D Gauss-Fit

    # return list of x,y positions of local maxima
    return np.asarray(np.where(maxima)).T; 


# --- self-test -------------------------------------------------------------
if __name__ == '__main__':
  import tifffile as tiff

  # read test image
  image   = tiff.imread("tests/graphene_flower_filtered.tif");
  FH      = FindCenters(image,verbosity=3);

  plt.show();
