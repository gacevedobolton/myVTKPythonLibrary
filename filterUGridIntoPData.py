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

def filterUGridIntoPData(
        ugrid,
        only_trianlges=False,
        verbose=1):

    myVTK.myPrint(verbose, "*** filterUGridIntoPData ***")

    filter_geometry = vtk.vtkGeometryFilter()
    filter_geometry.SetInputData(ugrid)
    filter_geometry.Update()
    pdata = filter_geometry.GetOutput()

    if (only_trianlges):
        filter_triangle = vtk.vtkTriangleFilter()
        filter_triangle.SetInputData(pdata)
        filter_triangle.Update()
        pdata = filter_triangle.GetOutput()

    return pdata
