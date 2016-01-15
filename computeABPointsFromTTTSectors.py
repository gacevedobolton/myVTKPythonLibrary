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

def computeABPointsFromTTTSectors(
        ugrid_sectors,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeABPointsFromTTTSectors ***")

    n_points = ugrid_sectors.GetNumberOfPoints()
    n_cells = ugrid_sectors.GetNumberOfCells()

    n_csects = 12
    n_rsects = 3
    n_slices = n_points / (n_rsects+1) / (n_csects+1)
    myVTK.myPrint(verbose, "n_slices =", n_slices)

    zmin = ugrid_sectors.GetPoint(0)[2]
    zmax = ugrid_sectors.GetPoint(ugrid_sectors.GetNumberOfPoints()-1)[2]

    dist_btw_slices = abs(zmin-zmax)/(n_slices-1)
    myVTK.myPrint(verbose, "dist_btw_slices =", dist_btw_slices)

    A = ugrid_sectors.GetPoints().GetPoint(0)
    B = ugrid_sectors.GetPoints().GetPoint(6)
    C = ugrid_sectors.GetPoints().GetPoint(3)
    D = ugrid_sectors.GetPoints().GetPoint(9)

    #print A
    #print B
    #print C
    #print D

    Px = ((A[0]*B[1]-A[1]*B[0])*(C[0]-D[0])-(A[0]-B[0])*(C[0]*D[1]-C[1]*D[0]))/((A[0]-B[0])*(C[1]-D[1])-(A[1]-B[1])*(C[0]-D[0]))
    Py = ((A[0]*B[1]-A[1]*B[0])*(C[1]-D[1])-(A[1]-B[1])*(C[0]*D[1]-C[1]*D[0]))/((A[0]-B[0])*(C[1]-D[1])-(A[1]-B[1])*(C[0]-D[0]))

    #print Px
    #print Py

    A = [Px, Py, zmin]
    B = [Px, Py, zmax]

    myVTK.myPrint(verbose, "A =", A)
    myVTK.myPrint(verbose, "B =", B)

    points_AB = vtk.vtkPoints()
    points_AB.InsertNextPoint(A)
    points_AB.InsertNextPoint(B)

    cells_AB = vtk.vtkCellArray()
    cell_AB  = vtk.vtkVertex()
    cell_AB.GetPointIds().SetId(0, 0)
    cells_AB.InsertNextCell(cell_AB)
    cell_AB.GetPointIds().SetId(0, 1)
    cells_AB.InsertNextCell(cell_AB)

    AB_ugrid = vtk.vtkUnstructuredGrid()
    AB_ugrid.SetPoints(points_AB)
    AB_ugrid.SetCells(vtk.VTK_VERTEX, cells_AB)

    return points_AB



