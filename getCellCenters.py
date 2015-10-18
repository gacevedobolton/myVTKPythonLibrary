#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def getCellCenters(
        mesh,
        verbose=1):

    myVTK.myPrint(verbose, "*** getCellCenters ***")

    filter_cell_centers = vtk.vtkCellCenters()
    filter_cell_centers.SetInputData(mesh)
    filter_cell_centers.Update()

    return filter_cell_centers.GetOutput()
