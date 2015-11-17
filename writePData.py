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

import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def writePData(
        pdata,
        filename,
        verbose=1):

    myVTK.myPrint(verbose, "*** writePData: " + filename + " ***")

    if ('vtk' in filename):
        pdata_writer = vtk.vtkPolyDataWriter()
    elif ('vtp' in filename):
        pdata_writer = vtk.vtkXMLPolyDataWriter()
    else:
        assert 0, "File must be .vtk or .vtp. Aborting."

    pdata_writer.SetFileName(filename)
    pdata_writer.SetInputData(pdata)
    pdata_writer.Update()
    pdata_writer.Write()
