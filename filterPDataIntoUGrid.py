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

def filterPDataIntoUGrid(
        pdata,
        verbose=1):

    myVTK.myPrint(verbose, "*** filterPDataIntoUGrid ***")

    filter_append = vtk.vtkAppendFilter()
    if (vtk.vtkVersion.GetVTKMajorVersion() >= 6):
        filter_append.SetInputData(pdata)
    else:
        filter_append.SetInput(pdata)
    filter_append.Update()

    return filter_append.GetOutput()
