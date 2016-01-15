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

def thresholdUGrid(
        ugrid_mesh,
        field_support,
        field_name,
        threshold_value,
        threshold_by_upper_or_lower,
        verbose=1):

    myVTK.myPrint(verbose, "*** thresholdUGrid ***")

    threshold = vtk.vtkThreshold()
    threshold.SetInputData(ugrid_mesh)
    if (field_support == "points"):
        association = vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS
    elif (field_support == "cells"):
        association = vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS
    threshold.SetInputArrayToProcess(0, 0, 0, association, field_name)
    if (threshold_by_upper_or_lower == "upper"):
        threshold.ThresholdByUpper(threshold_value)
    elif (threshold_by_upper_or_lower == "lower"):
        threshold.ThresholdByLower(threshold_value)
    threshold.Update()
    ugrid_thresholded_mesh = threshold.GetOutput()
    return ugrid_thresholded_mesh

def thresholdPData(pdata_mesh,
                   field_support,
                   field_name,
                   threshold_value,
                   threshold_by_upper_or_lower,
                   verbose=1):
    myVTK.myPrint(verbose, "*** thresholdPData ***")

    ugrid_thresholded_mesh = thresholdUGrid(pdata_mesh,
                                            field_support,
                                            field_name,
                                            threshold_value,
                                            threshold_by_upper_or_lower,
                                            False)
    pdata_thresholded_mesh = filterUGridIntoPData(ugrid_thresholded_mesh, False)
    return pdata_thresholded_mesh