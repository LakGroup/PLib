# vorntest.py
# script to perform voronoi analysis on Shreyasi's data

import numpy as np
from scipy import spatial as sp
import sys
from binLoad import binLoad
from vispy import gloo, app
import copy
from collections import defaultdict
import itertools

#FileLoc = 'E:/Data/rdemo/storm1.bin'
FileLoc = 'E:/Data/clusterT/raw/storm2_2527-1842_ch1.bin'
OutDict = binLoad(FileLoc)

# 2D-polyArea function
def PolyArea(X):
    return (0.5*np.abs(np.dot(X[:,0],np.roll(X[:,1],1))-
            np.dot(X[:,1],np.roll(X[:,0],1))))

# Pull out X Y coordinates only
Pos = OutDict['Localizations'][:,0:2]
Pos[:,1] = -Pos[:,1] # flip the y-axis
uPos = np.unique(Pos,axis=0) # make sure points are unique
nPoints = uPos.shape[0]

## start Voronoi analysis here
# only automated analysis for now
# perform monte carlo step here
# start with pixel measurements for now, convert at the end!
tri = sp.Delaunay(uPos)
vor = sp.Voronoi(uPos,qhull_options='Qbb Qc Qx')
Area = np.zeros(nPoints)
Vertices = dict()
for counter, index in enumerate(vor.point_region):
    vertIndex = vor.regions[index]
    Vertices[counter] = vor.vertices[vertIndex]
    Area[counter] = PolyArea(Vertices[counter])

# find the neighbors of the surrounding points
neiList=defaultdict(set)
for point in tri.vertices:
    for i,j in itertools.combinations(point,2):
        neiList[i].add(j)
        neiList[j].add(i)

AreaNeighbors = np.zeros(nPoints)
for ii in range(0,nPoints):
    AreaNeighbors[ii] = np.sum(Area[list(neiList[ii])])

from scipy.spatial import voronoi_plot_2d
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
#voronoi_plot_2d(vor)

#relArea = Area/max(Area)
# colorize
#for region, value in zip(vor.regions, relArea):
#    if not -1 in region:
#        polygon = [vor.vertices[i] for i in region]
#        plt.fill(*zip(*polygon), color='r',  alpha=value)

#for i,p in enumerate(uPos):
#    plt.text(p[0], p[1], '#%d' % i, ha='center')
#plt.show()

# lets plot this stuff

# find min/max values for normalization
minima = min(np.log(Area))
#maxima = max(Area)
maxima = max(np.log(Area))

# normalize chosen colormap
norm = mpl.colors.Normalize(vmin=minima, vmax=maxima, clip=True)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.seismic)

# plot Voronoi diagram, and fill finite regions with color mapped from speed value
voronoi_plot_2d(vor, show_points=False, show_vertices=False, s=1)
for r in range(len(vor.point_region)):
    region = vor.regions[vor.point_region[r]]
    if not -1 in region:
        polygon = [vor.vertices[i] for i in region]
        plt.fill(*zip(*polygon), color=mapper.to_rgba(np.log(Area[r])))
plt.show()
