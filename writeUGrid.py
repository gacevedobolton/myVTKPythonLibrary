#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def writeUGrid(
        ugrid,
        filename,
        verbose=1):

    myVTK.myPrint(verbose, "*** writeUGrid: " + filename + " ***")

    if ('vtk' in filename):
        ugrid_writer = vtk.vtkUnstructuredGridWriter()
    elif ('vtu' in filename):
        ugrid_writer = vtk.vtkXMLUnstructuredGridWriter()
    else:
        assert 0, "File must be .vtk or .vtu. Aborting."

    ugrid_writer.SetFileName(filename)
    ugrid_writer.SetInputData(ugrid)
    ugrid_writer.Update()
    ugrid_writer.Write()
