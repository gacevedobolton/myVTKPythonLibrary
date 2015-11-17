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

import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def getPointLocator(
        mesh,
        verbose=1):

    myVTK.myPrint(verbose, "*** getPointLocator ***")

    point_locator = vtk.vtkPointLocator()
    point_locator.SetDataSet(mesh)
    point_locator.Update()

    return point_locator
