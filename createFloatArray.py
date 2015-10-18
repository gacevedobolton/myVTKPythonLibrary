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

def createFloatArray(
        name,
        n_components=1,
        n_tuples=0,
        verbose=1):

    farray = vtk.vtkFloatArray()
    farray.SetName(name)
    farray.SetNumberOfComponents(n_components)
    farray.SetNumberOfTuples(n_tuples)

    return farray
