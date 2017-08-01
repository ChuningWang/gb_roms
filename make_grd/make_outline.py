'''
http://blog.thehumangeo.com/2014/05/12/drawing-boundaries-in-python/
'''

from shapely.ops import cascaded_union, polygonize
from shapely.geometry import Point
import shapely.geometry as geometry
from scipy.spatial import Delaunay
from matplotlib import path
import numpy as np
import math
import netCDF4 as nc

# --------------------------------------------------------------
def alpha_shape(coords, alpha):
    """
    Compute the alpha shape (concave hull) of a set
    of points.
    @param coords: Iterable container of coords.
    @param alpha: alpha value to influence the
        gooeyness of the border. Smaller numbers
        don't fall inward as much as larger numbers.
        Too large, and you lose everything!
    """
    if len(coords) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return geometry.MultiPoint(list(points)).convex_hull
    def add_edge(edges, edge_points, coords, i, j):
        """
        Add a line between the i-th and j-th points,
        if not in the list already
        """
        if (i, j) in edges or (j, i) in edges:
            # already added
            return
        edges.add( (i, j) )
        edge_points.append(coords[ [i, j] ])
    tri = Delaunay(coords)
    edges = set()
    edge_points = []
    # loop over triangles:
    # ia, ib, ic = indices of corner points of the
    # triangle
    for ia, ib, ic in tri.vertices:
        pa = coords[ia]
        pb = coords[ib]
        pc = coords[ic]
        # Lengths of sides of triangle
        a = math.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
        b = math.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
        c = math.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
        # Semiperimeter of triangle
        s = (a + b + c)/2.0
        # Area of triangle by Heron's formula
        area = math.sqrt(s*(s-a)*(s-b)*(s-c))
        circum_r = a*b*c/(4.0*area)
        # Here's the radius filter.
        #print circum_r
        if circum_r < 1.0/alpha:
            add_edge(edges, edge_points, coords, ia, ib)
            add_edge(edges, edge_points, coords, ib, ic)
            add_edge(edges, edge_points, coords, ic, ia)
    m = geometry.MultiLineString(edge_points)
    triangles = list(polygonize(m))
    return cascaded_union(triangles), edge_points

def plot_polygon(polygon):
    from descartes import PolygonPatch
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111)
    margin = .3
    x_min, y_min, x_max, y_max = polygon.bounds
    ax.set_xlim([x_min-margin, x_max+margin])
    ax.set_ylim([y_min-margin, y_max+margin])
    patch = PolygonPatch(polygon, fc='#999999',
                                          ec='#000000', fill=True,
                                          zorder=-1)
    ax.add_patch(patch)
    return fig

# --------------------------------------------------------------
import read_host_info
sv = read_host_info.read_host_info()
bathy_dir = sv['in_dir']

fh = nc.Dataset(bathy_dir + 'bathy_usgs.nc', 'r')
lon1 = fh.variables['lon'][:][2::5, 2::5]
lat1 = fh.variables['lat'][:][2::5, 2::5]
z1 = fh.variables['z'][:][2::5, 2::5]
fh.close()

msk = ~z1.mask
lon1 = lon1[msk]
lat1 = lat1[msk]
z1 = z1[msk]

# # step 1, draw boundary for usgs grid, only need once
# ct = len(z1)
# coords = np.array([(lon1[i], lat1[i]) for i in range(ct)])
# 
# concave_hull, edge_points = alpha_shape(coords, alpha=150.)
# bdry = np.array(concave_hull.boundary.coords[:])
# 
# # save boundary into text
# np.savetxt('bdry_usgs.txt', bdry)

bdry = np.loadtxt('bdry_usgs.txt')

# step 2, load noaa data, pick points out of usgs grid.
fh = nc.Dataset(bathy_dir + 'bathy_noaa.nc', 'r')
lon2 = fh.variables['lon'][:]
lat2 = fh.variables['lat'][:]
z2 = fh.variables['z'][:]
fh.close()

p0 = [(bdry[i, 0], bdry[i, 1]) for i in range(len(bdry[:, 0]))]

p = path.Path(p0)
pc = ~p.contains_points(np.array([lon2, lat2]).T) 

lon2 = lon2[pc]
lat2 = lat2[pc]
z2 = z2[pc]

lon = np.concatenate((lon1, lon2))
lat = np.concatenate((lat1, lat2))
z = np.concatenate((z1, z2))

# --------------------------------------------------------------
# lon = lon1[:8485]
# lat = lat1[:8485]
# z = z1[:8485]

ct = len(z)
coords = np.array([(lon[i], lat[i]) for i in range(ct)])

# --------------------------------------------------------------
concave_hull, edge_points = alpha_shape(coords, alpha=125.)
bdry = np.array(concave_hull.boundary.coords[:])

np.savetxt('bdry.txt', bdry)

# # --------------------------------------------------------------
# fh = nc.Dataset('bdry.nc', 'w')
# fh.createDimension('pts')
# fh.createVariable('lon', 'd', ('pts'))
# fh.createVariable('lat', 'd', ('pts'))
# fh.variables['lon'][:] = bdry[:, 0]
# fh.variables['lat'][:] = bdry[:, 1]
# fh.close()
