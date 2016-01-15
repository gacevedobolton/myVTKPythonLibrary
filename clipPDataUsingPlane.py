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

import sys
import vtk

import myVTKPythonLibrary as myVTK
from mat_vec_tools import *

########################################################################

def clipPDataUsingPlane(
        pdata_mesh,
        plane_C,
        plane_N,
        verbose=1):

    myVTK.myPrint(verbose, "*** clipPDataUsingPlane ***")

    plane = vtk.vtkPlane()
    plane.SetOrigin(plane_C)
    plane.SetNormal(plane_N)

    clip = vtk.vtkClipPolyData()
    clip.SetClipFunction(plane)
    clip.GenerateClippedOutputOn()
    clip.SetInputData(pdata_mesh)
    clip.Update()
    clipped0 = clip.GetOutput(0)
    clipped1 = clip.GetOutput(1)

    if (clipped0.GetNumberOfPoints() > clipped1.GetNumberOfPoints()):
        return clipped0, clipped1
    else:
        return clipped1, clipped0
