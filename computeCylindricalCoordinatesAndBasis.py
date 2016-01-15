#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2016                               ###
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

def computeCylindricalCoordinatesAndBasis(
        points,
        points_AB,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeCylindricalCoordinatesAndBasis ***")

    assert (points_AB.GetNumberOfPoints() >= 2), "points_AB must have at least two points. Aborting."
    point_A = numpy.array([0.]*3)
    point_B = numpy.array([0.]*3)
    points_AB.GetPoint(                              0, point_A)
    points_AB.GetPoint(points_AB.GetNumberOfPoints()-1, point_B)
    eL  = point_B - point_A
    eL /= numpy.linalg.norm(eL)
    if (verbose >= 2): print "eL =", eL

    point_C = point_A+numpy.array([1.,0.,0.])
    #point_C = numpy.array(points.GetPoint(0))

    AC  = point_C - point_A
    AD  = numpy.cross(eL, AC)
    AD /= numpy.linalg.norm(AD)
    AC  = numpy.cross(AD, eL)

    n_points = points.GetNumberOfPoints()

    farray_r = myVTK.createFloatArray("r", 1, n_points)
    farray_c = myVTK.createFloatArray("c", 1, n_points)
    farray_l = myVTK.createFloatArray("l", 1, n_points)

    farray_eR = myVTK.createFloatArray("eR", 3, n_points)
    farray_eC = myVTK.createFloatArray("eC", 3, n_points)
    farray_eL = myVTK.createFloatArray("eL", 3, n_points)

    for k_point in xrange(n_points):
        if (verbose >= 2): print "k_point =", k_point

        point = points.GetPoint(k_point)

        if (verbose >= 2): print "point =", point

        eR  = point - point_A
        eC  = numpy.cross(eL, eR)
        eC /= numpy.linalg.norm(eC)
        eR  = numpy.cross(eC, eL)

        farray_eR.InsertTuple(k_point, eR)
        farray_eC.InsertTuple(k_point, eC)
        farray_eL.InsertTuple(k_point, eL)

        r = numpy.dot(point - point_A, eR)
        farray_r.InsertTuple(k_point, [r])

        c  = math.atan2(numpy.dot(eR, AD), numpy.dot(eR, AC))
        c += (c<0.)*(2*math.pi)
        farray_c.InsertTuple(k_point, [c])

        l = numpy.dot(point - point_A, eL)
        farray_l.InsertTuple(k_point, [l])

    return (farray_r,
            farray_c,
            farray_l,
            farray_eR,
            farray_eC,
            farray_eL)

########################################################################

def addCylindricalCoordinatesAndBasis(
        ugrid,
        points_AB,
        verbose=1):

    myVTK.myPrint(verbose, "*** addCylindricalCoordinatesAndBasis ***")

    points = ugrid.GetPoints()
    (farray_r,
     farray_c,
     farray_l,
     farray_eR,
     farray_eC,
     farray_eL) = computeCylindricalCoordinatesAndBasis(
        points=points,
        points_AB=points_AB,
        verbose=verbose-1)

    ugrid.GetPointData().AddArray(farray_r)
    ugrid.GetPointData().AddArray(farray_c)
    ugrid.GetPointData().AddArray(farray_l)
    ugrid.GetPointData().AddArray(farray_eR)
    ugrid.GetPointData().AddArray(farray_eC)
    ugrid.GetPointData().AddArray(farray_eL)

    cell_centers = myVTK.getCellCenters(
        mesh=ugrid,
        verbose=verbose-1)
    (farray_r,
     farray_c,
     farray_l,
     farray_eR,
     farray_eC,
     farray_eL) = computeCylindricalCoordinatesAndBasis(
        points=cell_centers,
        points_AB=points_AB,
        verbose=verbose-1)

    ugrid.GetCellData().AddArray(farray_r)
    ugrid.GetCellData().AddArray(farray_c)
    ugrid.GetCellData().AddArray(farray_l)
    ugrid.GetCellData().AddArray(farray_eR)
    ugrid.GetCellData().AddArray(farray_eC)
    ugrid.GetCellData().AddArray(farray_eL)
