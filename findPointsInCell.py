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

def findPointsInCell(points,
                     cell,
                     verbose=1):

    ugrid_cell = vtk.vtkUnstructuredGrid()
    ugrid_cell.SetPoints(cell.GetPoints())
    cell = vtk.vtkHexahedron()
    for k_point in xrange(8): cell.GetPointIds().SetId(k_point, k_point)
    cell_array_cell = vtk.vtkCellArray()
    cell_array_cell.InsertNextCell(cell)
    ugrid_cell.SetCells(vtk.VTK_HEXAHEDRON, cell_array_cell)

    geometry_filter = vtk.vtkGeometryFilter()
    geometry_filter.SetInputData(ugrid_cell)
    geometry_filter.Update()
    cell_boundary = geometry_filter.GetOutput()

    pdata_points = vtk.vtkPolyData()
    pdata_points.SetPoints(points)

    enclosed_points_filter = vtk.vtkSelectEnclosedPoints()
    enclosed_points_filter.SetSurfaceData(cell_boundary)
    enclosed_points_filter.SetInputData(pdata_points)
    enclosed_points_filter.Update()

    points_in_cell = [k_point for k_point in xrange(points.GetNumberOfPoints()) if enclosed_points_filter.GetOutput().GetPointData().GetArray('SelectedPoints').GetTuple(k_point)[0]]
    return points_in_cell
