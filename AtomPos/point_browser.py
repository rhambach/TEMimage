"""
  Simple GUI for manually specifying points on a 2D image

  USAGE
    An example can be found at the end of this file and can be executed
    using 'python point_browser.py'

  Copyright (c) 2013, rhambach. 
    This file is part of the FitAtomPos package and released
    under the MIT-Licence. See LICENCE file for details.
"""
import numpy as np
import matplotlib.pylab as plt
from matplotlib.widgets import Button, RadioButtons

class ImageBrowser(object):
  """
  ImageBrowser shows a grayscale image in a separate window
  """

  def __init__(self,image,imginfo={},fnext=None,verbosity=0):
    """
    image  ... 2D Array for background image
    imginfo... (opt) dictionary with parameters of the image:
               'desc'     ... image description (default: '')
               'filename' ... filename of the image (default: '')
               'xperchan' ... scaling for x-axis (default: 1)
               'yperchan' ... scaling for y-axis (default: 1)
               'xunits'   ... unit for x-axis (default: 'px')
               'yunits'   ... unit for y-axis (default: 'px')
               'atom'     ... atom species (default: 'C') 
    fnext  ... (opt) function(self,event) connected to the Next-button
    verbosity. (opt) quiet (0), verbose (3), debug (4)
    """

    self.image  = np.asarray(image);
    self.verbosity=verbosity;

    # set default values for imginfo
    self.imginfo= {'desc': "", 'filename':"", 'atom':'C', \
                   'xperchan':1., 'yperchan':1., 'xunits':'px', 'yunits':'px'};
    self.imginfo.update(imginfo);
    Imin = self._px2ic(-0.5,-0.5);
    Imax = self._px2ic(self.image.shape[0]-0.5, self.image.shape[1]-0.5); 
    self.imginfo['extent'] = [Imin[0],Imax[0], Imin[1],Imax[1]];
        
    # open new figure and draw image
    fig = plt.figure(); self.fig = fig;
    self.axis = self.fig.add_subplot(111); 
    self.fig.subplots_adjust(right=0.85);
    self.AxesImage=None;
    self.axis.set_title('ImageBrowser: %s' % self.imginfo['desc']);
    self._update_image(); 
   
    # add buttons
    axStyle = fig.add_axes([0.85, 0.1, 0.1, 0.08]);
    self.rbStyle = RadioButtons(axStyle,["blur","raw"]);
    self.rbStyle.on_clicked(self.ChangeStyle);
    if fnext is not None:
      axNext  = self.fig.add_axes([0.85, 0.2, 0.1, 0.04]);    
      self.bNext   = Button(axNext,'Next');
      self.bNext.on_clicked(lambda *args: fnext(self,*args));


  def ChangeStyle(self,label):
    """ 
    Change style of background image 
    """
    if label=="raw":
      self.AxesImage.set_interpolation("nearest");
    elif label=="blur":
      self.AxesImage.set_interpolation("bilinear");
    self._update();


  def _update_image(self):
    " redraw image "
    if self.AxesImage is not None: self.AxesImage.remove();
    self.AxesImage = self.axis.imshow(self.image.T,cmap=plt.gray(),\
        origin='lower',aspect='equal',extent=self.imginfo['extent']);
    self.axis.set_xlabel("x [%s]" % self.imginfo['xunits']);
    self.axis.set_ylabel("y [%s]" % self.imginfo['yunits']);
    self.axis.set_xlim(*self.imginfo['extent'][0:2]);
    self.axis.set_ylim(*self.imginfo['extent'][2:4]);
    self._update();

  def _update(self):
    plt.draw();

  def _ic2px(self,x,y):
    "convert image coordinates to pixel positions"
    return (x/self.imginfo['xperchan'], y/self.imginfo['yperchan']);

  def _px2ic(self,x,y):
    "convert pixel positions to image coordinates"
    return (x*self.imginfo['xperchan'], y*self.imginfo['yperchan']);





class PointBrowser(ImageBrowser):
  """
  PointBroswer allows to manually add or remove points in a 2D image
  """

  def __init__(self,image,points,imginfo={},fnext=None,verbosity=0):
    """
    image  ... 2D Array for background image
    points ... initial list of peak positions shape(N,2)
    imginfo... (opt) dictionary with image parameters (see ImageBrowser)
    self.bDual.on_clicked(self.Dual);
    verbosity. (opt) quiet (0), verbose (3), debug (4)
    """

    # init ImageBrowser
    super(PointBrowser,self).__init__(image,imginfo,fnext,verbosity);
    self.axis.set_title('PointBrowser: %s' % self.imginfo['desc']);

    # draw point list
    self.points = np.asarray(points);
    assert( self.points.ndim==2 and self.points.shape[1]==2 );
    self.Line2D, = self.axis.plot(self.points[:,0], self.points[:,1],\
                    'ro',picker=5);
    #self.axis.set_xlim(*self.imginfo['extent'][0:2]);
    #self.axis.set_ylim(*self.imginfo['extent'][2:4]);

    # add buttons
    axEdit  = self.fig.add_axes([0.85, 0.85,0.1, 0.04]);
    axSave  = self.fig.add_axes([0.85, 0.8, 0.1, 0.04]);
    axLoad  = self.fig.add_axes([0.85, 0.75,0.1, 0.04]);
    self.bEdit   = Button(axEdit,'Edit');
    self.bSave   = Button(axSave,'Save');
    self.bLoad   = Button(axLoad,'Load');
    self.bEdit.on_clicked(self.ToggleEdit);
    self.bSave.on_clicked(self.Save);
    self.bLoad.on_clicked(self.Load);

    self._update();
    

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
    self._update();


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

    self._update_points();

    
  def _update_points(self):
    self.Line2D.set_data(self.points.T);
    self._update();

  def Save(self,event):
    import tkFileDialog
    filename = tkFileDialog.asksaveasfilename(
         filetypes=(('XYZ File', '*.xyz'),
                    ('Text File','*.txt')));
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
        OUT.write("Point positions in (%s,%s) for image '%s'\n" \
            %(self.imginfo['xunits'],self.imginfo['yunits'],    \
              self.imginfo['filename']));
        # data for xyz file (in image coordinates)
        for x,y in self.points:
          OUT.write("%2s  %8f  %8f  0.000000\n" % (self.imginfo['atom'], x, y));
        OUT.close();
    except:
      from tkMessageBox import showerror
      showerror("Save error", "Could not write point list to file '%s'" % filename);

  def Load(self,event):
    """
    supported file formats: 
      *.txt:  list of x,y positions in pixels
      *.xyz:  absolute coordinates
    """
    import tkFileDialog
    filename = tkFileDialog.askopenfilename(
        filetypes=(('XYZ File', '*.xyz'),
                    ('Text File','*.txt')));
    if not filename: return;
    ext = filename.split('.')[-1];
    try:
      if self.verbosity>1: print("LOADING point list from file '%s'" % filename);    
      if ext=='txt':
        self.points=np.genfromtxt(filename,dtype=float);
      elif ext=='xyz':
        self.points=np.genfromtxt(filename,usecols=(1,2),dtype=float,skip_header=2);
    except:
      from tkMessageBox import showerror
      showerror("Load error", "could not read point list from file '%s'" % filename);
    self._update_points();





class TilingBrowser(ImageBrowser):
  """
  TilingBroswer allows to manually add or remove lines in a 2D image
  """

  def __init__(self,image,tiling,imginfo={},fnext=None,verbosity=0):
    """
    image  ... 2D Array for background image
    tiling ... Tiling instance
    imginfo... (opt) dictionary with image parameters (see ImageBrowser)
    self.bDual.on_clicked(self.Dual);
    verbosity. (opt) quiet (0), verbose (3), debug (4)
    """
    import matplotlib.collections as c;

    # init ImageBrowser
    super(TilingBrowser,self).__init__(image,imginfo,fnext,verbosity);
    self.axis.set_title('TilingBrowser: %s' % self.imginfo['desc']);
    
    # draw edges
    self.tiling= tiling;
    segments = [self.tiling.points[e] for e in self.tiling.edges];
    self.LineCol = c.LineCollection(segments);
    self.LineCol.set_color('blue');
    self.axis.add_collection(self.LineCol);
    
    # add buttons
    axEdit  = self.fig.add_axes([0.85, 0.85,0.1, 0.04]);
    axSave  = self.fig.add_axes([0.85, 0.8, 0.1, 0.04]);
    axLoad  = self.fig.add_axes([0.85, 0.75,0.1, 0.04]);
    self.bEdit   = Button(axEdit,'Edit');
    self.bSave   = Button(axSave,'Save');
    self.bLoad   = Button(axLoad,'Load');
    self.bEdit.on_clicked(self.ToggleEdit);
    self.bSave.on_clicked(self.Save);
    self.bLoad.on_clicked(self.Load);

    self._update();
    

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
    self._update();


  def __onpick(self,event):
    """
    implements the following actions upon a mouse click event
    1. Left mouse click:  flip edge
    """

    x = event.mouseevent.xdata
    y = event.mouseevent.ydata
    if self.verbosity>3: print x,y,event.artist;

    # 1. Left mouse click: flip corresponding edge
    if event.mouseevent.button==1:
      if event.artist <> self.axis:   return True;
      inside,prop = self.LineCol.contains(event.mouseevent);
      if inside: # flip at least one edge
        if self.verbosity>2:
          print "FLIP edge", prop['ind'][0], " at: ", x,y
        self.tiling.flip(prop['ind'][0]);
      self._update_tiling();

    
  def _update_tiling(self):
    segments = [self.tiling.points[e] for e in self.tiling.edges];
    self.LineCol.set_segments(segments);
    self._update();

  def Save(self,event):
    import tkFileDialog
    filename = tkFileDialog.asksaveasfilename(
         filetypes=(('Dump File','*.pkl'),));
    if not filename: return;
    ext = filename.split('.')[-1];
    try:
      if self.verbosity>1: print("SAVING edges to file '%s'" % filename)
      if ext=='pkl':
        import pickle
        OUT = open(filename, 'wb');
        pickle.dump(self.tiling, OUT);
        OUT.close();
    except:
      from tkMessageBox import showerror
      showerror("Save error", "Could not write list of edges to file '%s'" % filename);

  def Load(self,event):
    """
    supported file formats: 
      *.txt:  list of x,y positions in pixels
      *.xyz:  absolute coordinates
    """
    import tkFileDialog
    filename = tkFileDialog.askopenfilename(
         filetypes=(('Dump File','*.pkl'),));
    if not filename: return;
    ext = filename.split('.')[-1];
    try:
      if self.verbosity>1: print("LOADING edges from file '%s'" % filename);    
      if ext=='pkl':
        import pickle;
        IN = open(filename, 'rb');
        self.tiling=pickle.load(IN);
        IN.close();
    except:
      from tkMessageBox import showerror
      showerror("Load error", "could not read point list from file '%s'" % filename);
    self._update_tiling();



# --- self-test -------------------------------------------------------------
if __name__ == '__main__':

  # ImageBrowser
  x,y = np.ogrid[0.:0.9:.02, 0.:0.8:.01];
  im  = np.exp(-((x-0.6)**2+(y-0.3)**2)/0.2);
  info= {'desc':'Test Image',
         'xperchan': 0.02, 'yperchan': 0.01,
         'xunits': 'nm', 'yunits': 'nm'};
  IB =ImageBrowser(im,info,verbosity=4);
  
  # PointBrowser
  pt  = np.round(np.random.rand(100, 2),decimals=2).tolist();
  pt.append([0.22,0.02]); 
  PB=PointBrowser(im,pt,info,verbosity=4);

  # TilingBrowser
  from tiling import *
  TB=TilingBrowser(im,Tiling(pt),info,verbosity=4);

  plt.show();

