#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import os
import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def readPData(
        filename,
        verbose=1):

    myVTK.myPrint(verbose, "*** readPData: " + filename + " ***")

    if ('vtk' in filename):
        pdata_reader = vtk.vtkPolyDataReader()
    elif ('vtp' in filename):
        pdata_reader = vtk.vtkXMLPolyDataReader()
    else:
        assert 0, "File must be .vtk or .vtp. Aborting."

    assert (os.path.isfile(filename)), "Wrong filename. Aborting."

    pdata_reader.SetFileName(filename)
    pdata_reader.Update()
    pdata = pdata_reader.GetOutput()

    if (verbose):
        print "n_points = " + str(pdata.GetNumberOfPoints())
        print "n_verts = " + str(pdata.GetNumberOfVerts())
        print "n_lines = " + str(pdata.GetNumberOfLines())
        print "n_polys = " + str(pdata.GetNumberOfPolys())
        print "n_strips = " + str(pdata.GetNumberOfStrips())

    return pdata
