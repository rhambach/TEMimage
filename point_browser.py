"""
     Simple GUI for manually specifying points on a 2D image

   VERSION
     $Id: dp_mdf.py 1604 2011-11-15 13:04:35Z hambach $

   USAGE
     An example can be found at the end of this file and can be executed
     using 'python point_browser.py'

   COPYRIGHT
     Copyright (c) 2011, Ralf Hambach. All rights reserved.
     Use of this source code is governed by a BSD-style license that 
     can be found in the LICENSE file.
"""
__version__ = "$Revision: 367 $"
__author__ = "Ralf Hambach"
__date__ = "$Date: 2011-11-15 13:58:13 +0100 (Tue, 15 Nov 2011) $"

import numpy as np
import matplotlib.pylab as plt
from matplotlib.widgets import Button, RadioButtons

class PointBrowser:
  """
  PointBroswer allows to manually add or remove points in a 2D image
    image  ... 2D Array for background image
    imginfo... dictionary with parameters of the image, required:
               'desc'     ... image description (default: '')
               'filename' ... filename of the image (default: '')
               'xperchan' ... scaling for x-axis (default: 1)
               'yperchan' ... scaling for y-axis (default: 1)
               'xunit'    ... unit for x-axis (default: 'px')
               'yunit'    ... unit for y-axis (default: 'px')
               'atom'     ... atom species (default: 'C')
    points ... initial list of peak positions shape(N,2)
    verbosity. (opt) quiet (0), verbose (3), debug (4)
  """

  def __init__(self,image,imginfo,points,verbosity=0):
    
    self.image  = np.asarray(image);
    # set default values for imginfo
    self.imginfo= {'desc': "", 'filename':"", 'atom':'C', \
                   'xperchan':1, 'yperchan':1, 'xunit':'px', 'yunit':'px'};
    self.imginfo.update(imginfo);
    self.points = np.asarray(points);
    assert( self.points.ndim==2 and self.points.shape[1]==2 );
    self.verbosity=verbosity;
        
    # open new figure with buttons
    fig = plt.figure(); self.fig = fig; 
    plt.title(self.imginfo['desc']);
    self.axis = fig.add_subplot(111); 
    self.fig.subplots_adjust(right=0.85);
    self.AxesImage = self.axis.imshow(self.image,cmap=plt.gray());
    self.Line2D,   = self.axis.plot(self.points[:,0], self.points[:,1],'ro',picker=5);
    self.axis.set_xlim(0,self.image.shape[0]-1);
    self.axis.set_ylim(self.image.shape[1]-1,0);

    # add buttons
    axStyle = fig.add_axes([0.85, 0.1, 0.1, 0.08]);
    axEdit  = fig.add_axes([0.85, 0.85,0.1, 0.04]);
    axSave  = fig.add_axes([0.85, 0.8, 0.1, 0.04]);
    axLoad  = fig.add_axes([0.85, 0.75,0.1, 0.04]);
    axDual  = fig.add_axes([0.85, 0.7, 0.1, 0.04]);
    self.rbStyle = RadioButtons(axStyle,["blur","raw"]);
    self.bEdit   = Button(axEdit,'Edit');
    self.bSave   = Button(axSave,'Save');
    self.bLoad   = Button(axLoad,'Load');
    self.bDual   = Button(axDual,'Dual');
    self.rbStyle.on_clicked(self.ChangeStyle);
    self.bEdit.on_clicked(self.ToggleEdit);
    self.bSave.on_clicked(self.Save);
    self.bLoad.on_clicked(self.Load);
    self.bDual.on_clicked(self.Dual);

    plt.show();
    
  def ChangeStyle(self,label):
    """ 
    Change style of background image 
    """
    if label=="raw":
      self.AxesImage.set_interpolation("nearest");
    elif label=="blur":
      self.AxesImage.set_interpolation("bilinear");
    plt.draw();

  def ToggleEdit(self,event):
    """ 
    Toggle interactive mode for adding/removing data points 
    """
    # Toggle ON
    if self.bEdit.label.get_text() == 'Edit':
      self.bEdit.label.set_text('Hold');
      self.axis.set_picker(True);
      self.cid_ToggleEdit = \
        self.fig.canvas.mpl_connect('pick_event', self.__onpick);
    # Toggle OFF
    else:
      self.bEdit.label.set_text('Edit');
      self.axis.set_picker(False);
      self.fig.canvas.mpl_disconnect(self.cid_ToggleEdit);
    plt.draw();

  def __onpick(self,event):
    """
    implements the following actions upon a mouse click event
    1. Right mouse click: remove clicked point
    2. Left mouse click:  add point at the given position
    """

    x = event.mouseevent.xdata
    y = event.mouseevent.ydata
    if self.verbosity>3: print x,y,event.artist;

    # 1. Right mouse click: remove point
    if event.mouseevent.button==3:
      if event.artist <> self.Line2D: return True;
      if self.verbosity>2: 
        print "REMOVE point ", event.ind, " at: ", self.points[event.ind]
      self.points=np.delete(self.points,event.ind,axis=0);
        
    # 2. Left mouse click: add point
    elif event.mouseevent.button==1:
      if event.artist <> self.axis:   return True;
      if self.verbosity>2:
        print "ADD point at: ", [x,y]
      self.points=np.append(self.points,[[x,y]],axis=0);

    self.__update();

    
  def __update(self):
    self.Line2D.set_data(self.points.T);
    plt.draw();

  def Save(self,event):
    import tkFileDialog
    filename = tkFileDialog.asksaveasfilename(
         filetypes=(('Text File','*.txt'),
                    ('XYZ File', '*.xyz')));
    if not filename: return;
    ext = filename.split('.')[-1];
    try:
      if self.verbosity>1: print("SAVING point list to file '%s'" % filename)
      if ext=='txt':
        # TODO: add header
        np.savetxt(filename,self.points);
      elif ext=='xyz':
        OUT=open(filename,'w');
        # header of xyz file
        OUT.write("%d\n"%len(self.points));
        OUT.write("Point positions in %s for image '%s'\n" \
            %(self.imginfo['xunit'],self.imginfo['filename']));
        # data for xyz file (in image coordinates)
        for x,y in self.points:
          x,y = self.__px2ic(x,y);
          OUT.write("%2s  %8f  %8f  0.000000\n" % (self.imginfo['atom'], x, y));
        OUT.close();
    except:
      from tkMessageBox import showerror
      showerror("Save error", "Could not write point list to file '%s'" % filename);
      raise


  def Load(self,event):
    """
    supported file formats: 
      *.txt:  list of x,y positions in pixels
      *.xyz:  absolute coordinates
    """
    import tkFileDialog
    filename = tkFileDialog.askopenfilename(
         filetypes=(('Text File','*.txt'),
                    ("XYZ file", "*.xyz")));
    if not filename: return;
    ext = filename.split('.')[-1];
    try:
      if self.verbosity>1: print("LOADING point list from file '%s'" % filename);    
      if ext=='txt':
        self.points=np.genfromtxt(filename,dtype=float);
      elif ext=='xyz':
        pts=np.genfromtxt(filename,usecols=(1,2),dtype=float,skip_header=2);
        self.points=np.asarray(self.__ic2px(pts[:,0],pts[:,1])).T;
    except:
      print("ERROR: could not read point list from file '%s'" % filename);
      raise
    self.__update();

  def __ic2px(self,x,y):
    "convert image coordinates to pixel positions"
    return (x/self.imginfo['xperchan'], y/self.imginfo['yperchan']);

  def __px2ic(self,x,y):
    "convert pixel positions to image coordinates"
    return (x*self.imginfo['xperchan'], y*self.imginfo['yperchan']);

  def Dual(self,event):
    """
    Raise new instance with dual points generated from the centers 
    of a Delaunay-triangulation
    """
    from scipy.spatial import Delaunay
    tri = Delaunay(self.points);
    dual = [];
    for edges in tri.vertices:
      dual.append(np.mean(self.points[edges],axis=0));
    info=self.imginfo.copy(); info['desc']+=", dual points";
    self.dualInstance = PointBrowser(self.image, info, dual);


# --- self-test -------------------------------------------------------------
if __name__ == '__main__':

  x,y = np.ogrid[0.:1.:.01, 0.:1.:.01];
  im  = np.exp(-((x-0.5)**2+(y-0.3)**2)/0.2);
  pt  = np.random.rand(100, 2)*100;
  info= {'desc':'Test Image','xperchan': 0.1};
  PB=PointBrowser(im,info,pt,verbosity=4);


