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

import os
import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def readSTL(
        filename,
        verbose=1):

    myVTK.myPrint(verbose, "*** readSTL: " + filename + " ***")

    assert (os.path.isfile(filename)), "Wrong filename. Aborting."

    stl_reader = vtk.vtkSTLReader()
    stl_reader.SetFileName(filename)
    stl_reader.Update()
    pdata_mesh = stl_reader.GetOutput()

    if (verbose):
        n_points = pdata_mesh.GetNumberOfPoints()
        print 'n_points =', n_points

        n_cells = pdata_mesh.GetNumberOfCells()
        print 'n_cells =', n_cells

    return pdata_mesh
