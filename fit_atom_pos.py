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
from hyperspy.drawing.image import plot_image_peaks
import point_browser as pb
from tiling import *

infile="tests/graphene_defect_filtered.dm3";
rawfile="tests/graphene_defect_raw.dm3";
IN = load(infile);
INRAW=load(rawfile);

# 1. Fourier Filtering for graphene (ring filter)
# Todo


# 2. Peak Fitting using Hyperspy

# normalize data to interval [0,1000] (for two_dim_findpeaks)
data=np.nan_to_num(IN.data);
data-=np.min(data);
data*=1000/np.max(data);
peaks=two_dim_findpeaks(data,peak_width=10);
center,height = peaks[:,0:2], peaks[:,2];

# 3. Start Point Browser
a=IN.axes_manager.axes; 
info = { 'desc':     IN.mapped_parameters.title,
         'filename': infile,
         'xperchan': a[0].scale,
         'yperchan': a[1].scale,
         'xunit'   : a[0].units,
         'yunit'   : a[1].units,
         'atoms'   : 'C' };
# convert peak positions from px -> unit
center *= [a[0].scale, a[1].scale];


# 3. Dual
def Dual(event):
  """
  Raise new instance with dual points generated from the centers 
  of a Delaunay-triangulation
  """
  Centers=Tiling(PBcenter.points);
  Atoms  =Centers.get_dual();
  info   =PBcenter.imginfo.copy(); 
  info['desc']=INRAW.mapped_parameters.title+", dual points";
  #PBatoms=pb.PointBrowser(IN.data,Atoms.points,info);
  PBimage=pb.ImageBrowser(INRAW.data,info);
  ax=Atoms.plot_edges(PBimage.axis,fc='darkgray');
  Atoms.plot_tiles(ax,nvertices=[5],fc='green',alpha=0.3);
  Atoms.plot_tiles(ax,nvertices=[7],fc='blue',alpha=0.2);
  Atoms.plot_tiles(ax,nvertices=[1,2,3,4,8,9,10],fc='red');

  plt.show();


PBcenter=pb.PointBrowser(IN.data,center,info,fnext=Dual);
plt.show();
