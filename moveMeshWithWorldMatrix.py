#!/usr/bin/python
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

import argparse
import numpy
import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def moveMeshWithWorldMatrix(
        mesh,
        M,
        verbose=1):

    myVTK.myPrint(verbose, "*** moveMeshWithWorldMatrix ***")

    n_points = mesh.GetNumberOfPoints()

    P = numpy.array([0.]*4)

    for k_point in xrange(n_points):
        P[0:3] = mesh.GetPoints().GetPoint(k_point)
        P[3] = 1.
        #print P

        P = numpy.dot(M, P)
        #print new_P

        mesh.GetPoints().SetPoint(k_point, P[0:3])
