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
from mat_vec_tools import *

########################################################################

def maskImageWithMesh(
    image,
    mesh,
    filter_with_field=None,
    verbose=1):

    myVTK.myPrint(verbose, "*** maskImageWithMesh ***")

    n_points = image.GetNumberOfPoints()
    farray_scalars_image = image.GetPointData().GetArray("scalars") # note that the field is defined at the points, not the cells

    (cell_locator,
     closest_point,
     generic_cell,
     k_cell,
     subId,
     dist) = myVTK.getCellLocator(
        mesh=mesh,
        verbose=verbose-1)

    if (filter_with_field != None):
        field = mesh.GetCellData().GetArray(filter_with_field[0])
        field_values = filter_with_field[1]

    points = vtk.vtkPoints()

    farray_scalars = myVTK.createFloatArray(
        name="scalars",
        n_components=1)

    for k_point in xrange(n_points):
        point = image.GetPoint(k_point)

        k_cell = cell_locator.FindCell(point)
        if (k_cell == -1): continue
        if (filter_with_field != None) and (field.GetTuple(k_cell)[0] not in field_values): continue

        points.InsertNextPoint(point)
        farray_scalars.InsertNextTuple(farray_scalars_image.GetTuple(k_point))

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.GetPointData().AddArray(farray_scalars)
    myVTK.addVertices(
        ugrid)

    return ugrid











