"""
     Simple GUI for visualising and analysing FFT patterns

   USAGE
     An example can be found at the end of this file and can be executed
     using 'python fft_plot.py'

   COPYRIGHT
     Copyright (c) 2011, Ralf Hambach. All rights reserved.
     Use of this source code is governed by a BSD-style license that 
     can be found in the LICENSE file.
"""
__version__ = "$Revision: 15 $"
__author__ = "Ralf Hambach"
__date__ = "$Date: 2011-12-13 23:20:48 +0100 (Tue, 13 Dec 2011) $"

import numpy as np
import matplotlib.pylab as plt
from matplotlib.widgets import Button, RadioButtons, Slider
import wqplot.wq_stack  as wq

class FFTBrowser(wq.WQBrowser):
  """
  WQStackBrowser shows a series of grayscale images in a separate window
  """

  def __init__(self,image,imginfo={},**kwargs):
    """
    image  ... 2D Array 
    imginfo... (opt) dictionary with parameters of the image:
               'desc'     ... image description (default: '')
               'filename' ... filename of the image (default: '')
               'xlabel'   ... name for x-axis    (default: x)
               'ylabel'   ... name for y-axis    (default: y)
               'xperchan' ... scaling for x-axis (default: 1)
               'yperchan' ... scaling for y-axis (default: 1)
               'xunits'   ... unit for x-axis (default: 'px')
               'yunits'   ... unit for y-axis (default: 'px')
    verbosity. (opt) quiet (0), verbose (3), debug (4)
    futher options are passed to the imshow method of matplotlib
    """
    # init ImageBrowser
    super(FFTBrowser,self).__init__(image,imginfo,**kwargs);
    self.fig.subplots_adjust(right=0.85, bottom=0.2);

    # add Switches and Buttons
    axBasis = self.fig.add_axes([0.85, 0.5, 0.12, 0.05]);    
    self.bBasis = Button(axBasis,'Set Basis');
    self.bBasis.on_clicked(self._toggle_set_basis);
    self.bBasis.color_activated='red';
    self.BasisSelector=None;

    # finally draw image 
    self._reset_image();     

  def _toggle_set_basis(self,event):
    if self.BasisSelector is None:
      O = np.asarray(self.image.shape)/2;
      self.BasisSelector=BasisSelector(self,O=O);
    else:
      print self.BasisSelector;
      self.BasisSelector=None;
    # swap colors
    b=self.bBasis; 
    b.color,b.color_activated=b.color_activated,b.color   
    b.hovercolor=b.color;   

class BasisSelector():
  """
  BasisSelector interactively chooses a grid on top of a figure
  """
  def __init__(self,ImageBrowser,a1=None,a2=None,O=None):
    """
    
    """
    self.IB = ImageBrowser;
    self.a1 = np.asarray([0,10] if a1 is None else a1,dtype=float);
    self.a2 = np.asarray([10,0] if a2 is None else a2,dtype=float);
    self.O  = np.asarray([0,0.] if O  is None else O ,dtype=float);
    self.p0 = self.p = None;

    # draw grid
    self.points = np.asarray([(i*self.a1+j*self.a2+self.O) \
          for i in range(-5,6) for j in range(-5,6)]);
    self.Line2D,=self.IB.axis.plot(self.points[:,0],self.points[:,1],'or');
    #print 'init', self.Line2D, self.IB.axis.lines;
    #print 'test', self.Line2D in self.IB.axis.lines;

    # connect all events
    connect = self.IB.fig.canvas.mpl_connect;
    self.cidpress  = connect('button_press_event',  self.on_press)
    self.cidrelease= connect('button_release_event',self.on_release)
    self.cidmotion = connect('motion_notify_event', self.on_motion)

  def __del__(self):
    # disconnect all events
    self.IB.axis.lines.remove(self.Line2D);
    disconnect = self.IB.fig.canvas.mpl_disconnect;
    disconnect(self.cidpress);
    disconnect(self.cidrelease);
    disconnect(self.cidmotion);

  def __str__(self):
    return "Basis vectors: a1=(%g,%g)\n"%tuple(self.a1) \
          +"               a2=(%g,%g)  "%tuple(self.a2);

  def on_press(self,event):
    if event.inaxes != self.IB.axis: return
    self.p = self.p0 = event.xdata, event.ydata; 
                       # p0 is reference point for mouse motion
    #print 'press'

  def on_release(self,event):
    # transform using the last stored p and p0
    self.a1, self.a2 = self._transform_points([self.a1+self.O,self.a2+self.O])-self.O;
    print "basis vectors: a1: ",self.a1, ", a2: ", self.a2
    self.points = np.asarray([(i*self.a1+j*self.a2+self.O) \
          for i in range(-5,6) for j in range(-5,6)]);
    self.IB.axis.lines.remove(self.Line2D);
    self.Line2D,=self.IB.axis.plot(self.points[:,0],self.points[:,1],'or');

    self.p0=None;
    #print 'release'
 
  def on_motion(self,event):
    if self.p0 is None: return
    if event.inaxes != self.IB.axis: return
    self.p = event.xdata,event.ydata;  # change of mouse position from p0
    
    # calculate new points (rotation around self.O)
    points=self._transform_points(self.points);         
    self.Line2D.set_data(points.T);
    self.IB._update();

  def _transform_points(self,points):
    """ 
    apply a rotation/zoom to given points according to change
    of the vector p from p0
    """
    p0=self.p0-self.O; p=self.p-self.O;
    alpha=np.arctan2(*p) - np.arctan2(*p0);
    s = np.sin(alpha); c = np.cos(alpha);
    rot = np.asarray([[ c, -s], [ s, c]]);      # rotation matrix
    scale=np.linalg.norm(p)/np.linalg.norm(p0); # scaling  factor
    return scale*np.dot(points-self.O,rot)+self.O;



# -- main ----------------------------------------
if __name__ == '__main__':
  import tifffile as tiff

  filename="tests/TiO2_0eV.tif";
  img =tiff.imread(filename);
  diff= np.fft.fftshift(np.fft.fft2(img));
  info={'xlabel': 'y','ylabel': 'x'};
  fig3=FFTBrowser(np.abs(diff),info,interpolation='none');
  plt.show();
