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

def computeHelixTransverseSheetAngles(
        farray_eRR,
        farray_eCC,
        farray_eLL,
        farray_eF,
        farray_eS,
        farray_eN,
        use_new_definition=False,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeHelixTransverseSheetAngles ***")

    n_tuples = farray_eRR.GetNumberOfTuples()
    assert (farray_eCC.GetNumberOfTuples() == n_tuples)
    assert (farray_eLL.GetNumberOfTuples() == n_tuples)
    assert (farray_eF.GetNumberOfTuples() == n_tuples)
    assert (farray_eS.GetNumberOfTuples() == n_tuples)
    assert (farray_eN.GetNumberOfTuples() == n_tuples)

    farray_angle_helix = myVTK.createFloatArray("angle_helix", 1, n_tuples)
    farray_angle_trans = myVTK.createFloatArray("angle_trans", 1, n_tuples)
    farray_angle_sheet = myVTK.createFloatArray("angle_sheet", 1, n_tuples)

    for k_tuple in xrange(n_tuples):
        eRR = numpy.array(farray_eRR.GetTuple(k_tuple))
        eCC = numpy.array(farray_eCC.GetTuple(k_tuple))
        eLL = numpy.array(farray_eLL.GetTuple(k_tuple))

        eF  = numpy.array(farray_eF.GetTuple(k_tuple))
        eF -= numpy.dot(eF, eRR) * eRR
        eF /= numpy.linalg.norm(eF)
        angle_helix = math.copysign(1., numpy.dot(eF, eCC)) * math.asin(min(1., max(-1., numpy.dot(eF, eLL)))) * (180./math.pi)
        farray_angle_helix.InsertTuple(k_tuple, [angle_helix])

        eF  = numpy.array(farray_eF.GetTuple(k_tuple))
        eF -= numpy.dot(eF, eLL) * eLL
        eF /= numpy.linalg.norm(eF)
        angle_trans = math.copysign(-1., numpy.dot(eF, eCC)) * math.asin(min(1., max(-1., numpy.dot(eF, eRR)))) * (180./math.pi)
        farray_angle_trans.InsertTuple(k_tuple, [angle_trans])

        #if (use_new_definition):
            #assert 0, "TODO"
        #else:
            #assert 0, "TODO"

    return (farray_angle_helix,
            farray_angle_trans,
            farray_angle_sheet)

########################################################################

def addHelixTransverseSheetAngles(
        ugrid,
        type_of_support="cell",
        use_new_definition=False,
        verbose=1):

    myVTK.myPrint(verbose, "*** addHelixTransverseSheetAngles ***")

    if (type_of_support == "cell"):
        ugrid_data = ugrid.GetCellData()
    elif (type_of_support == "point"):
        ugrid_data = ugrid.GetPointData()

    farray_eRR = ugrid_data.GetArray("eRR")
    farray_eCC = ugrid_data.GetArray("eCC")
    farray_eLL = ugrid_data.GetArray("eLL")

    farray_eF = ugrid_data.GetArray("eF")
    farray_eS = ugrid_data.GetArray("eS")
    farray_eN = ugrid_data.GetArray("eN")

    (farray_angle_helix,
     farray_angle_trans,
     farray_angle_sheet) = computeHelixTransverseSheetAngles(
        farray_eRR=farray_eRR,
        farray_eCC=farray_eCC,
        farray_eLL=farray_eLL,
        farray_eF=farray_eF,
        farray_eS=farray_eS,
        farray_eN=farray_eN,
        use_new_definition=use_new_definition,
        verbose=verbose-1)

    ugrid_data.AddArray(farray_angle_helix)
    ugrid_data.AddArray(farray_angle_trans)
    ugrid_data.AddArray(farray_angle_sheet)

    return (farray_angle_helix,
            farray_angle_trans,
            farray_angle_sheet)
