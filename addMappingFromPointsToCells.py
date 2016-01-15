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
import math
import numpy

import myVTKPythonLibrary as myVTK

########################################################################

def addMappingFromPointsToCells(
        ugrid_points,
        ugrid_cells,
        verbose=1):

    myVTK.myPrint(verbose, "*** addMappingFromPointsToCells ***")

    n_points = ugrid_points.GetNumberOfPoints()
    n_cells = ugrid_cells.GetNumberOfCells()
    #print "n_points = " + str(n_points)
    #print "n_cells = " + str(n_cells)

    (cell_locator,
     closest_point,
     generic_cell,
     k_cell,
     subId,
     dist) = getCellLocator(
         mesh=ugrid_cells,
         verbose=verbose-1)

    iarray_k_cell = createIntArray("k_cell", 1, n_points)

    for k_point in xrange(n_points):
        point = ugrid_points.GetPoint(k_point)

        cell_locator.FindClosestPoint(
            point,
            closest_point,
            generic_cell,
            k_cell,
            subId,
            dist)
        #k_cell = cell_locator.FindCell(point)

        iarray_k_cell.InsertTuple(k_point, [k_cell])
        #print "k_point = " + str(k_point)
        #print "k_cell = " + str(k_cell)

    #ugrid_points.GetPointData().AddArray(iarray_k_cell)
    ugrid_points.GetCellData().AddArray(iarray_k_cell)
