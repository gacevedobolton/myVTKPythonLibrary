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

def addPDataNormals(
        pdata,
        verbose=1):

    myVTK.myPrint(verbose, "*** addPDataNormals ***")

    poly_data_normals = vtk.vtkPolyDataNormals()
    poly_data_normals.ComputePointNormalsOff()
    poly_data_normals.ComputeCellNormalsOn()
    poly_data_normals.SetInputData(pdata)
    poly_data_normals.Update()
    pdata = poly_data_normals.GetOutput()

    mesh_center = numpy.array(pdata.GetCenter())
    #print mesh_center

    cell_centers = myVTK.getCellCenters(
        mesh=pdata,
        verbose=verbose-1)

    cnt_pos = 0
    cnt_neg = 0
    for k_cell in xrange(pdata.GetNumberOfCells()):
        cell_center = cell_centers.GetPoint(k_cell)
        outward  = cell_center-mesh_center
        outward /= numpy.linalg.norm(outward)
        normal = pdata.GetCellData().GetNormals().GetTuple(k_cell)
        proj = numpy.dot(outward, normal)
        if (proj > 0): cnt_pos += 1
        else:          cnt_neg += 1
    #print cnt_pos
    #print cnt_neg

    if (cnt_neg > cnt_pos):
        poly_data_normals.FlipNormalsOn()
        poly_data_normals.Update()
        pdata = poly_data_normals.GetOutput()

    return pdata
