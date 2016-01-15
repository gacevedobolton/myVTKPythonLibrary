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

import os
import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def readUGrid(
        filename,
        verbose=1):

    myVTK.myPrint(verbose, "*** readUGrid: " + filename + " ***")

    if ('vtk' in filename):
        ugrid_reader = vtk.vtkUnstructuredGridReader()
    elif ('vtu' in filename):
        ugrid_reader = vtk.vtkXMLUnstructuredGridReader()
    else:
        assert 0, "File must be .vtk or .vtu. Aborting."

    assert (os.path.isfile(filename)), "Wrong filename. Aborting."

    ugrid_reader.SetFileName(filename)
    ugrid_reader.Update()
    ugrid = ugrid_reader.GetOutput()

    if (verbose):
        n_points = ugrid.GetNumberOfPoints()
        print 'n_points =', n_points

        n_cells = ugrid.GetNumberOfCells()
        print 'n_cells =', n_cells

    return ugrid
