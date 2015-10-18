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

def writeSTL(
        pdata,
        filename,
        verbose=1):

    myVTK.myPrint(verbose, "*** writeSTL: " + filename + " ***")

    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetFileName(filename)
    stl_writer.SetInputData(pdata)
    stl_writer.Update()
    stl_writer.Write()
