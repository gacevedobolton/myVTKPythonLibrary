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

import math
import numpy
import vtk
import os
import sys

import myVTKPythonLibrary as myVTK

########################################################################

def computeRegionsForBiV(
        points,
        pdata_endLV,
        pdata_endRV,
        pdata_epi,
        verbose=0):

    myVTK.myPrint(verbose, "*** computeRegionsForBiV ***")

    myVTK.myPrint(verbose, "Initializing cell locators...")

    (cell_locator_endLV,
     closest_point_endLV,
     generic_cell,
     cellId_endLV,
     subId,
     dist_endLV) = myVTK.getCellLocator(
         mesh=pdata_endLV,
         verbose=verbose-1)
    (cell_locator_endRV,
     closest_point_endRV,
     generic_cell,
     cellId_endRV,
     subId,
     dist_endRV) = myVTK.getCellLocator(
         mesh=pdata_endRV,
         verbose=verbose-1)
    (cell_locator_epi,
     closest_point_epi,
     generic_cell,
     cellId_epi,
     subId,
     dist_epi) = myVTK.getCellLocator(
         mesh=pdata_epi,
         verbose=verbose-1)

    n_points = points.GetNumberOfPoints()

    iarray_region = myVTK.createIntArray("region_id", 1, n_points)

    for k_point in range(n_points):
        point = numpy.array(points.GetPoint(k_point))
        cell_locator_endLV.FindClosestPoint(
            point,
            closest_point_endLV,
            generic_cell,
            cellId_endLV,
            subId,
            dist_endLV)
        cell_locator_endRV.FindClosestPoint(
            point,
            closest_point_endRV,
            generic_cell,
            cellId_endRV,
            subId,
            dist_endRV)
        cell_locator_epi.FindClosestPoint(
            point,
            closest_point_epi,
            generic_cell,
            cellId_epi,
            subId,
            dist_epi)

        if   (dist_endRV == max(dist_endLV, dist_endRV, dist_epi)):
            iarray_region.SetTuple(k_point, [0])
        elif (dist_epi == max(dist_endLV, dist_endRV, dist_epi)):
            iarray_region.SetTuple(k_point, [1])
        elif (dist_endLV == max(dist_endLV, dist_endRV, dist_epi)):
            iarray_region.SetTuple(k_point, [2])

    return iarray_region

########################################################################

def addRegionsToBiV(
        ugrid_mesh,
        pdata_endLV,
        pdata_endRV,
        pdata_epi,
        verbose=0):

    myVTK.myPrint(verbose, "*** addRegionsToBiV ***")

    points = ugrid_mesh.GetPoints()
    iarray_region = computeRegionsForBiV(
        points=points,
        pdata_endLV=pdata_endLV,
        pdata_endRV=pdata_endRV,
        pdata_epi=pdata_epi,
        verbose=verbose-1)
    ugrid_mesh.GetPointData().AddArray(iarray_region)

    cell_centers = myVTK.getCellCenters(
        mesh=ugrid_mesh,
        verbose=verbose-1)
    iarray_region = computeRegionsForBiV(
        points=cell_centers,
        pdata_endLV=pdata_endLV,
        pdata_endRV=pdata_endRV,
        pdata_epi=pdata_epi,
        verbose=verbose-1)
    ugrid_mesh.GetCellData().AddArray(iarray_region)

########################################################################

