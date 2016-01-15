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

def computePseudoProlateSpheroidalCoordinatesAndBasisForLV(
        points,
        farray_c,
        farray_l,
        farray_eL,
        pdata_end,
        pdata_epi,
        iarray_part_id=None,
        verbose=1):

    myVTK.myPrint(verbose, "*** computePseudoProlateSpheroidalCoordinatesAndBasisForLV ***")

    myVTK.myPrint(verbose, "Computing surface cell normals...")

    pdata_end = myVTK.addPDataNormals(
        pdata=pdata_end,
        verbose=verbose-1)
    pdata_epi = myVTK.addPDataNormals(
        pdata=pdata_epi,
        verbose=verbose)

    myVTK.myPrint(verbose, "Initializing surface cell locators...")

    (cell_locator_end,
     closest_point_end,
     generic_cell,
     cellId_end,
     subId,
     dist_end) = myVTK.getCellLocator(
         mesh=pdata_end,
         verbose=verbose-1)
    (cell_locator_epi,
     closest_point_epi,
     generic_cell,
     cellId_epi,
     subId,
     dist_epi) = myVTK.getCellLocator(
         mesh=pdata_epi,
         verbose=verbose-1)

    myVTK.myPrint(verbose, "Computing local prolate spheroidal directions...")

    n_points = points.GetNumberOfPoints()

    farray_rr = myVTK.createFloatArray("rr", 1, n_points)
    farray_cc = myVTK.createFloatArray("cc", 1, n_points)
    farray_ll = myVTK.createFloatArray("ll", 1, n_points)

    farray_eRR = myVTK.createFloatArray("eRR", 3, n_points)
    farray_eCC = myVTK.createFloatArray("eCC", 3, n_points)
    farray_eLL = myVTK.createFloatArray("eLL", 3, n_points)

    if (n_points == 0):
        return (farray_rr,
                farray_cc,
                farray_ll,
                farray_eRR,
                farray_eCC,
                farray_eLL)

    c_lst = [farray_c.GetTuple(k_point)[0] for k_point in xrange(n_points)]
    c_min = min(c_lst)
    c_max = max(c_lst)

    l_lst = [farray_l.GetTuple(k_point)[0] for k_point in xrange(n_points)]
    l_min = min(l_lst)
    l_max = max(l_lst)

    for k_point in xrange(n_points):
        if (iarray_part_id is not None) and (int(iarray_part_id.GetTuple(k_point)[0]) > 0):
            rr = 0.
            cc = 0.
            ll = 0.
            eRR = [1.,0.,0.]
            eCC = [0.,1.,0.]
            eLL = [0.,0.,1.]

        else:
            point = numpy.array(points.GetPoint(k_point))
            cell_locator_end.FindClosestPoint(
                point,
                closest_point_end,
                generic_cell,
                cellId_end,
                subId,
                dist_end)
            cell_locator_epi.FindClosestPoint(
                point,
                closest_point_epi,
                generic_cell,
                cellId_epi,
                subId,
                dist_epi)

            rr = dist_end/(dist_end+dist_epi)

            c = farray_c.GetTuple(k_point)[0]
            cc = (c-c_min) / (c_max-c_min)

            l = farray_l.GetTuple(k_point)[0]
            ll = (l-l_min) / (l_max-l_min)

            normal_end = numpy.reshape(pdata_end.GetCellData().GetNormals().GetTuple(cellId_end), (3))
            normal_epi = numpy.reshape(pdata_epi.GetCellData().GetNormals().GetTuple(cellId_epi), (3))
            eRR  = (1.-rr) * normal_end + rr * normal_epi
            eRR /= numpy.linalg.norm(eRR)

            eL = numpy.reshape(farray_eL.GetTuple(k_point), (3))
            eCC  = numpy.cross(eL, eRR)
            eCC /= numpy.linalg.norm(eCC)

            eLL = numpy.cross(eRR, eCC)

        farray_rr.InsertTuple(k_point, [rr])
        farray_cc.InsertTuple(k_point, [cc])
        farray_ll.InsertTuple(k_point, [ll])
        farray_eRR.InsertTuple(k_point, eRR)
        farray_eCC.InsertTuple(k_point, eCC)
        farray_eLL.InsertTuple(k_point, eLL)

    return (farray_rr,
            farray_cc,
            farray_ll,
            farray_eRR,
            farray_eCC,
            farray_eLL)

########################################################################

def addPseudoProlateSpheroidalCoordinatesAndBasisToLV(
        ugrid,
        pdata_end,
        pdata_epi,
        verbose=1):

    myVTK.myPrint(verbose, "*** addPseudoProlateSpheroidalCoordinatesAndBasisToLV ***")

    (farray_rr,
     farray_cc,
     farray_ll,
     farray_eRR,
     farray_eCC,
     farray_eLL) = computePseudoProlateSpheroidalCoordinatesAndBasisForLV(
        points=ugrid.GetPoints(),
        farray_c=ugrid.GetPointData().GetArray("c"),
        farray_l=ugrid.GetPointData().GetArray("l"),
        farray_eL=ugrid.GetPointData().GetArray("eL"),
        pdata_end=pdata_end,
        pdata_epi=pdata_epi,
        iarray_part_id=ugrid.GetPointData().GetArray("part_id"),
        verbose=verbose-1)
    ugrid.GetPointData().AddArray(farray_rr)
    ugrid.GetPointData().AddArray(farray_cc)
    ugrid.GetPointData().AddArray(farray_ll)
    ugrid.GetPointData().AddArray(farray_eRR)
    ugrid.GetPointData().AddArray(farray_eCC)
    ugrid.GetPointData().AddArray(farray_eLL)

    cell_centers = myVTK.getCellCenters(
        mesh=ugrid,
        verbose=verbose-1)
    (farray_rr,
    farray_cc,
    farray_ll,
    farray_eRR,
    farray_eCC,
    farray_eLL) = computePseudoProlateSpheroidalCoordinatesAndBasisForLV(
        points=cell_centers,
        farray_c=ugrid.GetCellData().GetArray("c"),
        farray_l=ugrid.GetCellData().GetArray("l"),
        farray_eL=ugrid.GetCellData().GetArray("eL"),
        pdata_end=pdata_end,
        pdata_epi=pdata_epi,
        iarray_part_id=ugrid.GetCellData().GetArray("part_id"),
        verbose=verbose-1)
    ugrid.GetCellData().AddArray(farray_rr)
    ugrid.GetCellData().AddArray(farray_cc)
    ugrid.GetCellData().AddArray(farray_ll)
    ugrid.GetCellData().AddArray(farray_eRR)
    ugrid.GetCellData().AddArray(farray_eCC)
    ugrid.GetCellData().AddArray(farray_eLL)

########################################################################

def computePseudoProlateSpheroidalCoordinatesAndBasisForBiV(
        points,
        iarray_regions,
        farray_c,
        farray_l,
        farray_eL,
        pdata_endLV,
        pdata_endRV,
        pdata_epi,
        iarray_part_id=None,
        verbose=1):

    myVTK.myPrint(verbose, "*** computePseudoProlateSpheroidalCoordinatesAndBasisForBiV ***")

    myVTK.myPrint(verbose, "Computing surface cell normals...")

    pdata_endLV = myVTK.addPDataNormals(
        pdata=pdata_endLV,
        verbose=verbose-1)
    pdata_endRV = myVTK.addPDataNormals(
        pdata=pdata_endRV,
        verbose=verbose-1)
    pdata_epi = myVTK.addPDataNormals(
        pdata=pdata_epi,
        verbose=verbose)

    myVTK.myPrint(verbose, "Initializing surface cell locators...")

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

    myVTK.myPrint(verbose, "Computing local prolate spheroidal directions...")

    n_points = points.GetNumberOfPoints()

    farray_rr = myVTK.createFloatArray("rr", 1, n_points)
    farray_cc = myVTK.createFloatArray("cc", 1, n_points)
    farray_ll = myVTK.createFloatArray("ll", 1, n_points)

    farray_eRR = myVTK.createFloatArray("eRR", 3, n_points)
    farray_eCC = myVTK.createFloatArray("eCC", 3, n_points)
    farray_eLL = myVTK.createFloatArray("eLL", 3, n_points)

    c_lst_FWLV = numpy.array([farray_c.GetTuple(k_point)[0] for k_point in xrange(n_points) if (iarray_regions.GetTuple(k_point)[0] == 0)])
    (c_avg_FWLV, c_std_FWLV) = myVTK.computeMeanStddevAngles(
        angles=c_lst_FWLV,
        angles_in_degrees=False,
        angles_in_pm_pi=False)
    myVTK.myPrint(verbose, "c_avg_FWLV = " + str(c_avg_FWLV))
    c_lst_FWLV = (((c_lst_FWLV-c_avg_FWLV+math.pi)%(2*math.pi))-math.pi+c_avg_FWLV)
    c_min_FWLV = min(c_lst_FWLV)
    c_max_FWLV = max(c_lst_FWLV)
    myVTK.myPrint(verbose, "c_min_FWLV = " + str(c_min_FWLV))
    myVTK.myPrint(verbose, "c_max_FWLV = " + str(c_max_FWLV))

    c_lst_S = numpy.array([farray_c.GetTuple(k_point)[0] for k_point in xrange(n_points) if (iarray_regions.GetTuple(k_point)[0] == 1)])
    (c_avg_S, c_std_S) = myVTK.computeMeanStddevAngles(
        angles=c_lst_S,
        angles_in_degrees=False,
        angles_in_pm_pi=False)
    myVTK.myPrint(verbose, "c_avg_S = " + str(c_avg_S))
    c_lst_S = (((c_lst_S-c_avg_S+math.pi)%(2*math.pi))-math.pi+c_avg_S)
    c_min_S = min(c_lst_S)
    c_max_S = max(c_lst_S)
    myVTK.myPrint(verbose, "c_min_S = " + str(c_min_S))
    myVTK.myPrint(verbose, "c_max_S = " + str(c_max_S))

    c_lst_FWRV = numpy.array([farray_c.GetTuple(k_point)[0] for k_point in xrange(n_points) if (iarray_regions.GetTuple(k_point)[0] == 2)])
    (c_avg_FWRV, c_std_FWRV) = myVTK.computeMeanStddevAngles(
        angles=c_lst_FWRV,
        angles_in_degrees=False,
        angles_in_pm_pi=False)
    myVTK.myPrint(verbose, "c_avg_FWRV = " + str(c_avg_FWRV))
    c_lst_FWRV = (((c_lst_FWRV-c_avg_FWRV+math.pi)%(2*math.pi))-math.pi+c_avg_FWRV)
    c_min_FWRV = min(c_lst_FWRV)
    c_max_FWRV = max(c_lst_FWRV)
    myVTK.myPrint(verbose, "c_min_FWRV = " + str(c_min_FWRV))
    myVTK.myPrint(verbose, "c_max_FWRV = " + str(c_max_FWRV))

    l_lst = [farray_l.GetTuple(k_point)[0] for k_point in xrange(n_points)]
    l_min = min(l_lst)
    l_max = max(l_lst)

    for k_point in xrange(n_points):
        if (iarray_part_id is not None) and (int(iarray_part_id.GetTuple(k_point)[0]) > 0):
            rr = 0.
            cc = 0.
            ll = 0.
            eRR = [1.,0.,0.]
            eCC = [0.,1.,0.]
            eLL = [0.,0.,1.]

        else:
            point = numpy.array(points.GetPoint(k_point))
            region_id = iarray_regions.GetTuple(k_point)[0]

            if (region_id == 0):
                cell_locator_endLV.FindClosestPoint(
                    point,
                    closest_point_endLV,
                    generic_cell,
                    cellId_endLV,
                    subId,
                    dist_endLV)
                cell_locator_epi.FindClosestPoint(
                    point,
                    closest_point_epi,
                    generic_cell,
                    cellId_epi,
                    subId,
                    dist_epi)

                rr = dist_endLV/(dist_endLV+dist_epi)

                c = farray_c.GetTuple(k_point)[0]
                c = (((c-c_avg_FWLV+math.pi)%(2*math.pi))-math.pi+c_avg_FWLV)
                cc = (c-c_min_FWLV) / (c_max_FWLV-c_min_FWLV)

                l = farray_l.GetTuple(k_point)[0]
                ll = (l-l_min) / (l_max-l_min)

                normal_endLV = numpy.reshape(pdata_endLV.GetCellData().GetNormals().GetTuple(cellId_endLV), (3))
                normal_epi = numpy.reshape(pdata_epi.GetCellData().GetNormals().GetTuple(cellId_epi), (3))
                eRR  = (1.-rr) * normal_endLV + rr * normal_epi
                eRR /= numpy.linalg.norm(eRR)

                eL = numpy.reshape(farray_eL.GetTuple(k_point), (3))
                eCC  = numpy.cross(eL, eRR)
                eCC /= numpy.linalg.norm(eCC)

                eLL  = numpy.cross(eRR, eCC)
            elif (region_id == 1):
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

                rr = dist_endLV/(dist_endLV+dist_endRV)

                c = farray_c.GetTuple(k_point)[0]
                c = (((c-c_avg_S+math.pi)%(2*math.pi))-math.pi+c_avg_S)
                cc = (c-c_min_S) / (c_max_S-c_min_S)

                l = farray_l.GetTuple(k_point)[0]
                ll = (l-l_min) / (l_max-l_min)

                normal_endLV = numpy.reshape(pdata_endLV.GetCellData().GetNormals().GetTuple(cellId_endLV), (3))
                normal_endRV = numpy.reshape(pdata_endRV.GetCellData().GetNormals().GetTuple(cellId_endRV), (3))
                eRR  = (1.-rr) * normal_endLV - rr * normal_endRV
                eRR /= numpy.linalg.norm(eRR)

                eL = numpy.reshape(farray_eL.GetTuple(k_point), (3))
                eCC  = numpy.cross(eL, eRR)
                eCC /= numpy.linalg.norm(eCC)

                eLL = numpy.cross(eRR, eCC)
            if (region_id == 2):
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

                rr = dist_endRV/(dist_endRV+dist_epi)

                c = farray_c.GetTuple(k_point)[0]
                c = (((c-c_avg_FWRV+math.pi)%(2*math.pi))-math.pi+c_avg_FWRV)
                cc = (c-c_min_FWRV) / (c_max_FWRV-c_min_FWRV)

                l = farray_l.GetTuple(k_point)[0]
                ll = (l-l_min) / (l_max-l_min)

                normal_endRV = numpy.reshape(pdata_endRV.GetCellData().GetNormals().GetTuple(cellId_endRV), (3))
                normal_epi = numpy.reshape(pdata_epi.GetCellData().GetNormals().GetTuple(cellId_epi), (3))
                eRR  = (1.-rr) * normal_endRV + rr * normal_epi
                eRR /= numpy.linalg.norm(eRR)

                eL = numpy.reshape(farray_eL.GetTuple(k_point), (3))
                eCC  = numpy.cross(eL, eRR)
                eCC /= numpy.linalg.norm(eCC)

                eLL = numpy.cross(eRR, eCC)

        farray_rr.InsertTuple(k_point, [rr])
        farray_cc.InsertTuple(k_point, [cc])
        farray_ll.InsertTuple(k_point, [ll])
        farray_eRR.InsertTuple(k_point, eRR)
        farray_eCC.InsertTuple(k_point, eCC)
        farray_eLL.InsertTuple(k_point, eLL)

    return (farray_rr,
            farray_cc,
            farray_ll,
            farray_eRR,
            farray_eCC,
            farray_eLL)

########################################################################

def addPseudoProlateSpheroidalCoordinatesAndBasisToBiV(
        ugrid,
        pdata_endLV,
        pdata_endRV,
        pdata_epi,
        verbose=1):

    myVTK.myPrint(verbose, "*** addPseudoProlateSpheroidalCoordinatesAndBasisToBiV ***")

    (farray_rr,
     farray_cc,
     farray_ll,
     farray_eRR,
     farray_eCC,
     farray_eLL) = computePseudoProlateSpheroidalCoordinatesAndBasisForBiV(
        points=ugrid.GetPoints(),
        iarray_regions=ugrid.GetPointData().GetArray("region_id"),
        farray_c=ugrid.GetPointData().GetArray("c"),
        farray_l=ugrid.GetPointData().GetArray("l"),
        farray_eL=ugrid.GetPointData().GetArray("eL"),
        pdata_endLV=pdata_endLV,
        pdata_endRV=pdata_endRV,
        pdata_epi=pdata_epi,
        iarray_part_id=ugrid.GetPointData().GetArray("part_id"),
        verbose=verbose-1)
    ugrid.GetPointData().AddArray(farray_rr)
    ugrid.GetPointData().AddArray(farray_cc)
    ugrid.GetPointData().AddArray(farray_ll)
    ugrid.GetPointData().AddArray(farray_eRR)
    ugrid.GetPointData().AddArray(farray_eCC)
    ugrid.GetPointData().AddArray(farray_eLL)

    cell_centers = myVTK.getCellCenters(
        mesh=ugrid,
        verbose=verbose-1)
    (farray_rr,
     farray_cc,
     farray_ll,
     farray_eRR,
     farray_eCC,
     farray_eLL) = computePseudoProlateSpheroidalCoordinatesAndBasisForBiV(
        points=cell_centers,
        iarray_regions=ugrid.GetCellData().GetArray("region_id"),
        farray_c=ugrid.GetCellData().GetArray("c"),
        farray_l=ugrid.GetCellData().GetArray("l"),
        farray_eL=ugrid.GetCellData().GetArray("eL"),
        pdata_endLV=pdata_endLV,
        pdata_endRV=pdata_endRV,
        pdata_epi=pdata_epi,
        iarray_part_id=ugrid.GetCellData().GetArray("part_id"),
        verbose=verbose-1)
    ugrid.GetCellData().AddArray(farray_rr)
    ugrid.GetCellData().AddArray(farray_cc)
    ugrid.GetCellData().AddArray(farray_ll)
    ugrid.GetCellData().AddArray(farray_eRR)
    ugrid.GetCellData().AddArray(farray_eCC)
    ugrid.GetCellData().AddArray(farray_eLL)
