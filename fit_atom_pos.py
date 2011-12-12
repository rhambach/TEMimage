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

infile="tests/graphene_flower_filtered.dm3";
IN = load(infile);

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
pb.PointBrowser(IN.data,info,center);


# 3. Dual
def Dual(p):
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
