# fitting atomic positions in Fourier Filtered TEM Image
#
# 1. Fourier Filtering
# 2. Peak Fitting using HyperSpy
# 3. Delaunay Triangulation + Inversion to find Atom Positions

import numpy as np
import hyperspy;
from hyperspy.hspy          import load
from hyperspy.peak_char     import two_dim_findpeaks
from hyperspy.drawing.image import plot_image_peaks
import point_browser as pb

infile="tests/graphene_filter.dm3";
IN = load(infile);

# 1. Fourier Filtering for graphene (ring filter)
# Todo


# 2. Peak Fitting using Hyperspy
peaks=two_dim_findpeaks(IN.data);
center,height = peaks[:,0:2], peaks[:,2];
pb.PointBrowser(IN.data,center);
