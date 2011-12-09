"""
Calculates the dual of a Delaunay triangulation, namely a Voronoi diagram
naming convention:

Delaunay triangulation:
 vertex:   vertex
 triangle: face
 edge:     edge

Voronoi diagram:
 center:   vertex
 polygon:  face
 edge:     edge
"""

import numpy as np;
import matplotlib.pylab as plt;
from scipy.spatial import Delaunay;

class Voronoi:
  
  def __init__(self,tri):
    """
    calculate approximate Voronoi diagram from Delaunay triangulation tri
    """
    self.tri = tri;  # Delaunay triangulation

    # calculate the Dual
    self.points  = np.asarray([self.__get_point(t) for t in range(tri.nsimplex) ]);
    self.vertices= [self.__get_triangles(v) for v in range(tri.npoints) ];
    self.npoints = tri.nsimplex;
    self.nsimplex= tri.npoints;
  
  def __get_point(self,t):
    " center points of each triangle (dual to Delaunay face)"
    return np.mean(self.tri.points[self.tri.vertices[t]],axis=0)

  def __get_triangles(self,v):
    " return triangles that include Delaunay vertex v (sorted clockwise)"

    # search triangles that contain vertex v (t... triangle, n... which corner)
    t,n = np.where(tri.vertices == v);  
    # calculate vector from Delaunay vertex v to center of triangle t
    delta=self.points[t] - self.tri.points[v];
    # sort triagles using angle of vector to x-axis
    return t[ np.argsort(np.arctan2(*tuple(delta.T))) ];
  


# --- self-test -------------------------------------------------------------
if __name__ == '__main__':

  plt.gcf().clear();

  # Delaunay
  points = np.genfromtxt("test.txt");
  tri = Delaunay(points);
  #plt.plot(tri.points[:,0],tri.points[:,1],'bo');


  def plot_triangle(tri,t,color):
    vertices=tri.points[tri.vertices[t]].T.tolist();
    vertices[0].append(vertices[0][0]);
    vertices[1].append(vertices[1][0]);
    plt.plot(vertices[0],vertices[1],color);

  plot_triangle(tri,5,'r');
  for n in tri.neighbors[5]:
    if n<0: continue;
    plot_triangle(tri,n,'g');

  # Voronoi
  vor=Voronoi(tri);
  plt.plot(vor.points[:,0],vor.points[:,1],'r.');
  poly=vor.points[vor.vertices[100]];

  # Plot polygons
  N=20;
  
  polyhist=[ [[None,None]] for i in range(N) ];
  for poly in vor.vertices:
    if len(poly)>N-1: 
      print "WARNING: ignore polygon of size %d " % len(poly);
      continue;
    polyhist[len(poly)].extend(vor.points[poly].tolist()+[[None,None]]);

  # Irregular
  xirreg=[]; yirreg=[];
  for i in range(N):
    if i<5 or i>7:
      x,y = zip(*(polyhist[i]));
      xirreg.extend(x); yirreg.extend(y);
  plt.fill(xirreg,yirreg,facecolor='r',alpha=0.7, edgecolor='none')

  # Heptagons
  x,y = zip(*(polyhist[7]));
  plt.fill(x,y,facecolor='b',alpha=0.5, edgecolor='none')

  # Pentagons
  x,y = zip(*(polyhist[5]));
  plt.fill(x,y,facecolor='g',alpha=0.5, edgecolor='none')

  plt.show();

  

"""
This is maybe easiest to explain in code. The set of edges is:

	edges = []
	for i in xrange(x.nsimplex):
	    edges.append((x.vertices[i,0], x.vertices[i,1]))
	    edges.append((x.vertices[i,1], x.vertices[i,2]))
	    edges.append((x.vertices[i,2], x.vertices[i,0]))

This however counts each edge multiple times. To get around, that:

	edges = []
	for i in xrange(x.nsimplex):
	    if i > x.neighbors[i,2]:
	        edges.append((x.vertices[i,0], x.vertices[i,1]))
	    if i > x.neighbors[i,0]:
	        edges.append((x.vertices[i,1], x.vertices[i,2]))
	    if i > x.neighbors[i,1]:
	        edges.append((x.vertices[i,2], x.vertices[i,0]))
"""
