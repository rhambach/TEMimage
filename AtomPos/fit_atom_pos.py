# fitting atomic positions in Fourier Filtered TEM Image
#
# 1. Fourier Filtering
# 2. Peak Fitting using HyperSpy
# 3. Delaunay Triangulation 
# 4. Generate Atom Positions from Dual Graph

import numpy as np
import hyperspy;
from hyperspy.hspy          import load
from hyperspy.peak_char     import two_dim_findpeaks
import point_browser as pb
from tiling import *

def fit_peaks(data, xperchan=1, yperchan=1):
  " use hyperspy for finding centers of hexagons in image "

  # normalize data to [0,1000] (for two_dim_findpeaks)
  data=np.nan_to_num(data);
  data-=np.min(data);
  data*=1000/np.max(data);

  # peak fitting
  peaks=two_dim_findpeaks(data.T,peak_width=10);
  return peaks[:,0:2], peaks[:,2]; # centers, heights


def triangulation(PointBrowser, event):
  " raise LineBrowser to modify triangulation "
 
  # Delaunay-triangulation
  T=Tiling(PointBrowser.points); 
  # manually edit edges
  pb.TilingBrowser(INRAW.data,T,PointBrowser.imginfo,fnext=calc_dual);
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
    PBimage=pb.ImageBrowser(INRAW.data,info);
  elif event.button==2:
    PBimage=pb.PointBrowser(INRAW.data,dual.points,info,fnext=triangulation);
  if event.button==3:
    PBimage=pb.TilingBrowser(INRAW.data,dual,info,fnext=calc_dual);
  else:
    ax=dual.plot_edges(PBimage.axis,fc='darkgray');
    dual.plot_tiles(ax,nvertices=[5],fc='green',alpha=0.3);
    dual.plot_tiles(ax,nvertices=[7],fc='blue',alpha=0.2);
    dual.plot_tiles(ax,nvertices=[1,2,3,4,8,9,10],fc='red');

  plt.show();



# -- main ----------------------------------------


# read files (TODO: perform Fourier filter automatically)
#infile="tests/graphene_flower_filtered.dm3";
#rawfile="tests/graphene_flower_raw.dm3";
infile="tests/graphene_flower_filter.dm3";
rawfile="tests/graphene_flower_raw.dm3";
IN = load(infile);
INRAW=load(rawfile);

# image parameters
a=IN.axes_manager.axes; 
info = { 'desc':     IN.mapped_parameters.title,
         'filename': infile,
         'xperchan': a[0].scale,
         'yperchan': a[1].scale,
         'xunit'   : a[0].units,
         'yunit'   : a[1].units,
         'atoms'   : 'C' };

# fit + convert peak positions from px -> unit
centers, heights = fit_peaks(IN.data);
centers *= [info['xperchan'], info['yperchan']];

# manually edit center positions
PBcenter=pb.PointBrowser(IN.data,centers,info,fnext=triangulation);


plt.show();
