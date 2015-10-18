#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import math
import numpy
import vtk
import os
import sys

import myVTKPythonLibrary as myVTK

########################################################################

def computeSectorsForLV(
        farray_rr,
        farray_cc,
        farray_ll,
        n_r=1,
        n_c=1,
        n_l=1,
        iarray_part_id=None,
        verbose=0):

    myVTK.myPrint(verbose, "*** computeSectorsForLV ***")

    n_cells = farray_rr.GetNumberOfTuples()

    iarray_sector = myVTK.createIntArray("sector_id", 1, n_cells)

    for k_cell in range(n_cells):
        if (iarray_part_id is not None) and (int(iarray_part_id.GetTuple(k_cell)[0]) > 0):
            sector_id = -1

        else:
            rr = farray_rr.GetTuple(k_cell)[0]
            cc = farray_cc.GetTuple(k_cell)[0]
            ll = farray_ll.GetTuple(k_cell)[0]

            k_r = int(rr*n_r/1.000001)
            k_c = int(cc*n_c/1.000001)
            k_l = int((1.-ll)*n_l/1.000001)

            sector_id = k_l * n_c * n_r + k_c * n_r + k_r

        iarray_sector.SetTuple(k_cell, [sector_id])

    return iarray_sector

########################################################################

def addSectorsToLV(
        ugrid_mesh,
        n_r=1,
        n_c=1,
        n_l=1,
        verbose=0):

    myVTK.myPrint(verbose, "*** addSectorsToLV ***")

    iarray_sector = computeSectorsForLV(
        farray_rr=ugrid_mesh.GetCellData().GetArray("rr"),
        farray_cc=ugrid_mesh.GetCellData().GetArray("cc"),
        farray_ll=ugrid_mesh.GetCellData().GetArray("ll"),
        n_r=n_r,
        n_c=n_c,
        n_l=n_l,
        iarray_part_id=ugrid_mesh.GetCellData().GetArray("part_id"),
        verbose=verbose-1)

    ugrid_mesh.GetCellData().AddArray(iarray_sector)

########################################################################

def computeSectorsForBiV(
        iarray_regions,
        farray_rr,
        farray_cc,
        farray_ll,
        n_r=[1]*3,
        n_c=[1]*3,
        n_l=[1]*3,
        iarray_part_id=None,
        verbose=0):

    myVTK.myPrint(verbose, "*** computeSectorsForBiV ***")

    n_cells = iarray_regions.GetNumberOfTuples()
    assert (farray_rr.GetNumberOfTuples() == n_cells)
    assert (farray_cc.GetNumberOfTuples() == n_cells)
    assert (farray_ll.GetNumberOfTuples() == n_cells)

    iarray_sector = myVTK.createIntArray("sector_id", 1, n_cells)

    for k_cell in range(n_cells):
        if (iarray_part_id is not None) and (int(iarray_part_id.GetTuple(k_cell)[0]) > 0):
            sector_id = -1

        else:
            region_id = int(iarray_regions.GetTuple(k_cell)[0])

            rr = farray_rr.GetTuple(k_cell)[0]
            cc = farray_cc.GetTuple(k_cell)[0]
            ll = farray_ll.GetTuple(k_cell)[0]

            k_r = int(    rr *n_r[region_id]/1.000001)
            k_c = int(    cc *n_c[region_id]/1.000001)
            k_l = int((1.-ll)*n_l[region_id]/1.000001)

            sector_id = k_l * n_c[region_id] * n_r[region_id] \
                      + k_c * n_r[region_id] \
                      + k_r

            if (region_id >= 1):
                sector_id += n_r[0] * n_c[0] * n_l[0]
            if (region_id >= 2):
                sector_id += n_r[1] * n_c[1] * n_l[1]

        iarray_sector.SetTuple(k_cell, [sector_id])

    return iarray_sector

########################################################################

def addSectorsToBiV(
        ugrid_mesh,
        n_r=[1]*3,
        n_c=[1]*3,
        n_l=[1]*3,
        verbose=0):

    myVTK.myPrint(verbose, "*** addSectorsToBiV ***")

    iarray_sector = computeSectorsForBiV(
        iarray_regions=ugrid_mesh.GetCellData().GetArray("region_id"),
        farray_rr=ugrid_mesh.GetCellData().GetArray("rr"),
        farray_cc=ugrid_mesh.GetCellData().GetArray("cc"),
        farray_ll=ugrid_mesh.GetCellData().GetArray("ll"),
        n_r=n_r,
        n_c=n_c,
        n_l=n_l,
        iarray_part_id=ugrid_mesh.GetCellData().GetArray("part_id"),
        verbose=verbose-1)

    ugrid_mesh.GetCellData().AddArray(iarray_sector)

########################################################################
