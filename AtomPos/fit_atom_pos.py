# fitting atomic positions in Fourier Filtered TEM Image
#
# 1. Fourier Filtering
# 2. Peak Fitting
# 3. Delaunay Triangulation 
# 4. Generate Atom Positions from Dual Graph

import numpy as np
import tifffile as tiff
import point_browser as pb
from find_peaks import FindCenters
from tiling import *

def __add_info_text(ImageBrowser):
  " info about features of next-button"
  ax   = ImageBrowser.bNext.ax;
  text = "(use left/middle/right\n"+\
         "mouse button for\n"+\
         "preview/atoms/edges)";
  ax.text(-0.2, 1.5, text, fontsize=8, color='gray',
          transform=ax.transAxes);

def triangulation(PointBrowser, event):
  " raise LineBrowser to modify triangulation "
  if event.button<>1: return;
  # Delaunay-triangulation
  T=Tiling(PointBrowser.points); 
  # manually edit edges
  TB=pb.TilingBrowser(INRAW,T,PointBrowser.imginfo,fnext=calc_dual);
  __add_info_text(TB);
  
  plt.show();

def calc_dual(TilingBrowser,event):
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
    PBimage=pb.ImageBrowser(INRAW,info);
  elif event.button==2:
    PBimage=pb.PointBrowser(INRAW,dual.points,info,fnext=triangulation);
  if event.button==3:
    PBimage=pb.TilingBrowser(INRAW,dual,info,fnext=calc_dual);
    __add_info_text(PBimage);
  else:
    ax=dual.plot_edges(PBimage.axis,fc='darkgray');
    dual.plot_tiles(ax,nvertices=[5],fc='green',alpha=0.3);
    dual.plot_tiles(ax,nvertices=[7],fc='blue',alpha=0.2);
    dual.plot_tiles(ax,nvertices=[1,2,3,4,8,9,10],fc='red');

  plt.show();



# -- main ----------------------------------------


# read files (TODO: perform Fourier filter automatically)
infile="tests/graphene_flower_filtered.tif";
rawfile="tests/graphene_flower_raw.tif";
IN   =tiff.imread(infile);
INRAW=tiff.imread(rawfile);

# image parameters (TODO: read scale from Tiff?)
info = { 'desc':     infile.split('/')[-1],
         'filename': infile,
         'atoms'   : 'C' };

# fit + convert peak positions from px -> unit
FC = FindCenters(IN,imginfo=info,verbosity=3,fnext=triangulation);
plt.show();
