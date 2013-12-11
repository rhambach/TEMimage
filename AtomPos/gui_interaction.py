"""
  Organize GUI interaction for following steps:
    1. Fourier Filtering (ToDo)
    2. Peak Fitting
    3. Delaunay Triangulation 
    4. Generate Atom Positions from Dual Graph

  Copyright (c) 2013, rhambach. 
    This file is part of the TEMimage package and released
    under the MIT-Licence. See LICENCE file for details.
"""
import numpy as np
import tifffile as tiff
import point_browser as pb
from find_peaks import FindCenters
from tiling import *


class GUI(object):

  def __init__(self,infile,rawfile,verbosity=1):
    self.verbosity = verbosity;

    # read images
    self.IN   =tiff.imread(infile,verbosity=self.verbosity);
    self.INRAW=tiff.imread(rawfile,verbosity=self.verbosity);
  
    # image parameters (TODO: read scale from Tiff?)
    self.info = { 'desc':     infile.split('/')[-1],
                  'filename': infile,
                  'atoms'   : 'C' };

  def start(self): 
    FC = FindCenters(self.IN,imginfo=self.info,
          fnext=self.triangulation,verbosity=self.verbosity);
    plt.show(); # start GUI mainloop()

  def __add_info_text(self,ImageBrowser):
    " info about features of next-button"
    ax   = ImageBrowser.bNext.ax;
    text = "(use left/middle/right\n"+\
           "mouse button for\n"+\
           "preview/atoms/edges)";
    ax.text(-0.2, 1.5, text, fontsize=8, color='gray',
           transform=ax.transAxes);

  def triangulation(self, PointBrowser, event):
    " raise LineBrowser to modify triangulation "
    if event.button<>1: return;
    # Delaunay-triangulation
    T=Tiling(PointBrowser.points); 
    # manually edit edges
    TB=pb.TilingBrowser(self.INRAW,T,PointBrowser.imginfo,fnext=self.calc_dual);
    self.__add_info_text(TB);
    
    plt.show(); # raise figure window
  
  def calc_dual(self, TilingBrowser,event):
    """
    Raise new instance with dual points generated from the centers 
    of a Delaunay-triangulation
    """
    tri    = TilingBrowser.tiling;
    dual   = tri.get_dual();
    info   = TilingBrowser.imginfo.copy(); 
    info['desc']+=", dual points";
    # open different windows depending on mouse button
    if event.button==1:
      PBimage=pb.ImageBrowser(self.INRAW,info);
    elif event.button==2:
      PBimage=pb.PointBrowser(self.INRAW,dual.points,info,fnext=self.triangulation);
    if event.button==3:
      PBimage=pb.TilingBrowser(self.INRAW,dual,info,fnext=self.calc_dual);
      self.__add_info_text(PBimage);
    else:
      ax=dual.plot_edges(PBimage.axis,fc='darkgray');
      dual.plot_tiles(ax,nvertices=[5],fc='green',alpha=0.3);
      dual.plot_tiles(ax,nvertices=[7],fc='blue',alpha=0.2);
      dual.plot_tiles(ax,nvertices=[1,2,3,4,8,9,10],fc='red');
  
    plt.show(); # raise figure window
    


# -- main ----------------------------------------
if __name__=="__main__":

  # read files (TODO: perform Fourier filter automatically)
  infile="tests/graphene_flower_filtered.tif";
  rawfile="tests/graphene_flower_raw.tif";
  GUI(infile,rawfile).start();
   
