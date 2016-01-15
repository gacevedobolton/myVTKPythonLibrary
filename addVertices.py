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

import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def addVertices(
        ugrid,
        verbose=1):

    myVTK.myPrint(verbose, "*** addVertices ***")

    cell = vtk.vtkVertex()
    cell_array = vtk.vtkCellArray()

    n_points = ugrid.GetPoints().GetNumberOfPoints()
    for k_point in xrange(n_points):
        cell.GetPointIds().SetId(0, k_point)
        cell_array.InsertNextCell(cell)

    ugrid.SetCells(vtk.VTK_VERTEX, cell_array)

    n_arrays = ugrid.GetPointData().GetNumberOfArrays()
    for k_array in xrange(n_arrays):
        ugrid.GetCellData().AddArray(ugrid.GetPointData().GetArray(k_array))
