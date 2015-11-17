#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
### École Polytechnique, Palaiseau, France                           ###
###                                                                  ###
########################################################################

import math
import numpy

import myVTKPythonLibrary as myVTK

########################################################################

def computeCartesianCoordinates(
        points,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeCartesianCoordinates ***")

    n_points = points.GetNumberOfPoints()

    [xmin, xmax, ymin, ymax, zmin, zmax] = points.GetBounds()
    dx = xmax-xmin
    dy = ymax-ymin
    dz = zmax-zmin
    if (verbose >= 2): print "xmin = " + str(xmin)
    if (verbose >= 2): print "xmax = " + str(xmax)
    if (verbose >= 2): print "dx = " + str(dx)
    if (verbose >= 2): print "ymin = " + str(ymin)
    if (verbose >= 2): print "ymax = " + str(ymax)
    if (verbose >= 2): print "dy = " + str(dy)
    if (verbose >= 2): print "zmin = " + str(zmin)
    if (verbose >= 2): print "zmax = " + str(zmax)
    if (verbose >= 2): print "dz = " + str(dz)

    farray_norm_x = myVTK.createFloatArray("norm_x", 1, n_points)
    farray_norm_y = myVTK.createFloatArray("norm_y", 1, n_points)
    farray_norm_z = myVTK.createFloatArray("norm_z", 1, n_points)

    for k_point in xrange(n_points):
        if (verbose >= 2): print "k_point = " + str(k_point)

        point = points.GetPoint(k_point)
        if (verbose >= 2): print "point = " + str(point)

        norm_x = (point[0] - xmin) / dx
        norm_y = (point[1] - ymin) / dy
        norm_z = (point[2] - zmin) / dz

        farray_norm_x.InsertTuple(k_point, [norm_x])
        farray_norm_y.InsertTuple(k_point, [norm_y])
        farray_norm_z.InsertTuple(k_point, [norm_z])

    return (farray_norm_x,
            farray_norm_y,
            farray_norm_z)

########################################################################

def addCartesianCoordinates(
        ugrid,
        verbose=1):

    myVTK.myPrint(verbose, "*** addCartesianCoordinates ***")

    points = ugrid.GetPoints()
    (farray_norm_x,
     farray_norm_y,
     farray_norm_z) = computeCartesianCoordinates(
        points=points,
        verbose=verbose-1)

    ugrid.GetPointData().AddArray(farray_norm_x)
    ugrid.GetPointData().AddArray(farray_norm_y)
    ugrid.GetPointData().AddArray(farray_norm_z)

    cell_centers = myVTK.getCellCenters(
        mesh=ugrid,
        verbose=verbose-1)
    (farray_norm_x,
     farray_norm_y,
     farray_norm_z) = computeCartesianCoordinates(
        points=cell_centers,
        verbose=verbose-1)

    ugrid.GetCellData().AddArray(farray_norm_x)
    ugrid.GetCellData().AddArray(farray_norm_y)
    ugrid.GetCellData().AddArray(farray_norm_z)
