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

def mapDataSetToCellData(
        mesh_from,
        type_of_support,
        mesh_to,
        farray_names,
        type_of_mapping="PointsWithinRadius",
        radius_is_relative=True,
        radius=0.5,
        n_closest_points=3,
        threshold_dist=None,
        threshold_val_min=None,
        threshold_val_max=None,
        verbose=1):

    myVTK.myPrint(verbose, "*** mapDataSetToCellData ***")

    if (type_of_support == "point"):
        datapoints = mesh_from.GetPoints()
        dataset = mesh_from.GetPointData()
        point_locator = myVTK.getPointLocator(
            mesh_from,
            verbose=verbose-1)
    elif (type_of_support == "cell"):
        datapoints = myVTK.getCellCenters(
            mesh=mesh_from,
            verbose=verbose-1)
        dataset = mesh_from.GetCellData()
        point_locator = myVTK.getPointLocator(
            datapoints,
            verbose=verbose-1)
    else:
        assert (0)

    pdata_cell_centers_to = myVTK.getCellCenters(mesh_to)
    n_cells = mesh_to.GetNumberOfCells()

    farrays_avg = {}
    farrays_std = {}
    for farray_name in farray_names:
        assert (dataset.HasArray(farray_name)), "mesh has no array named " + farray_name + ". Aborting."

        farray_type = dataset.GetArray(farray_name).GetDataTypeAsString()
        farray_n_components = dataset.GetArray(farray_name).GetNumberOfComponents()
        farrays_avg[farray_name] = myVTK.createArray(farray_name+"_avg",
            farray_n_components,
            n_cells,
            farray_type)
        farrays_std[farray_name] = myVTK.createArray(farray_name+"_std",
            farray_n_components,
            n_cells,
            farray_type)

    points_within_radius = vtk.vtkIdList()

    for k_cell in xrange(n_cells):

        if (type_of_mapping == "ClosestPoints"):
            point_locator.FindClosestNPoints(
                n_closest_points,
                pdata_cell_centers_to.GetPoint(k_cell),
                points_within_radius)
        elif (type_of_mapping == "PointsWithinRadius"):
            if (radius_is_relative):
                l = (mesh_to.GetCell(k_cell).GetLength2())**(0.5)
                actual_radius = l*radius
            else:
                actual_radius = radius
            point_locator.FindPointsWithinRadius(
                actual_radius,
                pdata_cell_centers_to.GetPoint(k_cell),
                points_within_radius)
        else:
            assert (0)

        #points_within_radius = myVTK.findPointsInCell(mesh_from.GetPoints(), mesh_to.GetCell(k_cell))

        for farray_name in farray_names:
            if (points_within_radius.GetNumberOfIds()):
                values = [numpy.array(dataset.GetArray(farray_name).GetTuple(points_within_radius.GetId(k_id))) for k_id in xrange(points_within_radius.GetNumberOfIds()) if (threshold_dist == None) or (numpy.linalg.norm(numpy.array(datapoints.GetPoint(points_within_radius.GetId(k_id)))-numpy.array(pdata_cell_centers_to.GetPoint(k_cell))) < threshold_dist)]
                #print "values = " + str(values)
                if (threshold_val_min != None):
                    values = [value for value in values if (numpy.linalg.norm(value) > threshold_val_min)]
                if (threshold_val_max != None):
                    values = [value for value in values if (numpy.linalg.norm(value) < threshold_val_max)]
                #print "values = " + str(values)
                if (len(values)):
                    avg = numpy.mean(values, 0)
                    std = numpy.std(values, 0)
                else:
                    avg = [0]*farray_n_components
                    std = [0]*farray_n_components
            else:
                avg = [0]*farray_n_components
                std = [0]*farray_n_components
            farrays_avg[farray_name].InsertTuple(k_cell, avg)
            farrays_std[farray_name].InsertTuple(k_cell, std)

    for farray_name in farray_names:
        mesh_to.GetCellData().AddArray(farrays_avg[farray_name])
        mesh_to.GetCellData().AddArray(farrays_std[farray_name])
