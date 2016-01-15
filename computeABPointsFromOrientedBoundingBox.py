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

import numpy
import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def computeABPointsFromOrientedBoundingBox(
        mesh,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeABPointsFromOrientedBoundingBox ***")

    center = numpy.array(mesh.GetCenter())

    res_corner = numpy.array([0.]*3)
    res_max = numpy.array([0.]*3)
    res_mid = numpy.array([0.]*3)
    res_min = numpy.array([0.]*3)
    res_size = numpy.array([0.]*3)
    obb_tree = vtk.vtkOBBTree()
    obb_tree.ComputeOBB(
        mesh,
        res_corner,
        res_max,
        res_mid,
        res_min,
        res_size)
    if (verbose >= 2): print "res_corner =", res_corner
    if (verbose >= 2): print "res_max =", res_max
    if (verbose >= 2): print "res_mid =", res_mid
    if (verbose >= 2): print "res_min =", res_min
    if (verbose >= 2): print "res_size =", res_size

    point_A = res_corner + numpy.dot(center-res_corner, res_min/numpy.linalg.norm(res_min)) * res_min/numpy.linalg.norm(res_min) + numpy.dot(center-res_corner, res_mid/numpy.linalg.norm(res_mid)) * res_mid/numpy.linalg.norm(res_mid)
    point_B = point_A + res_max
    if (verbose >= 2): print "point_A =", point_A
    if (verbose >= 2): print "point_B =", point_B
    if (verbose >= 2): print "AB =", point_B-point_A

    points_AB = vtk.vtkPoints()
    points_AB.InsertNextPoint(point_A)
    points_AB.InsertNextPoint(point_B)
    #if (verbose >= 2): print points_AB

    return points_AB
