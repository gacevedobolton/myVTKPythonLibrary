#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import math
import numpy

import myVTKPythonLibrary as myVTK

########################################################################

def computeHelixAngles(
        farray_eRR,
        farray_eCC,
        farray_eLL,
        farray_eF,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeHelixAngles ***")

    n_tuples = farray_eRR.GetNumberOfTuples()
    assert (farray_eCC.GetNumberOfTuples() == n_tuples)
    assert (farray_eLL.GetNumberOfTuples() == n_tuples)
    assert (farray_eF.GetNumberOfTuples() == n_tuples)

    farray_angle_helix = myVTK.createFloatArray("angle_helix", 1, n_tuples)

    for k_tuple in xrange(n_tuples):
        eRR = numpy.array(farray_eRR.GetTuple(k_tuple))
        eCC = numpy.array(farray_eCC.GetTuple(k_tuple))
        eLL = numpy.array(farray_eLL.GetTuple(k_tuple))
        eF  = numpy.array(farray_eF.GetTuple(k_tuple))
        eF -= numpy.dot(eF, eRR) * eRR
        eF /= numpy.linalg.norm(eF)
        helix_angle = math.copysign(1., numpy.dot(eF, eCC)) * math.asin(min(1., max(-1., numpy.dot(eF, eLL)))) * (180./math.pi)
        farray_angle_helix.InsertTuple(k_tuple, [helix_angle])

    return farray_angle_helix

########################################################################

def addHelixAngles(
        ugrid,
        type_of_support="cell",
        verbose=1):

    myVTK.myPrint(verbose, "*** addHelixAngles ***")

    if (type_of_support == "cell"):
        ugrid_data = ugrid.GetCellData()
    elif (type_of_support == "point"):
        ugrid_data = ugrid.GetPointData()

    farray_eRR = ugrid_data.GetArray("eRR")
    farray_eCC = ugrid_data.GetArray("eCC")
    farray_eLL = ugrid_data.GetArray("eLL")

    farray_eF = ugrid_data.GetArray("eF")

    farray_angle_helix = computeHelixAngles(
        farray_eRR=farray_eRR,
        farray_eCC=farray_eCC,
        farray_eLL=farray_eLL,
        farray_eF=farray_eF,
        verbose=verbose-1)

    ugrid_data.AddArray(farray_angle_helix)

    return farray_angle_helix
