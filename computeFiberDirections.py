#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import sys
import math
import random
import numpy

import myVTKPythonLibrary as myVTK

########################################################################

def computeFiberDirections(
        farray_eRR,
        farray_eCC,
        farray_eLL,
        farray_angle_helix,
        angles_in_degrees=True,
        use_new_definition=False,
        shuffle_vectors=False,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeFiberDirections ***")

    n_tuples = farray_angle_helix.GetNumberOfTuples()

    farray_eF = myVTK.createFloatArray("eF", 3, n_tuples)
    farray_eS = myVTK.createFloatArray("eS", 3, n_tuples)
    farray_eN = myVTK.createFloatArray("eN", 3, n_tuples)

    for k_tuple in xrange(n_tuples):
        eRR = numpy.array(farray_eRR.GetTuple(k_tuple))
        eCC = numpy.array(farray_eCC.GetTuple(k_tuple))
        eLL = numpy.array(farray_eLL.GetTuple(k_tuple))

        assert (round(numpy.linalg.norm(eRR),1) == 1.0),\
            "|eRR| = " + str(numpy.linalg.norm(eRR)) + "≠ 1. Aborting"
        assert (round(numpy.linalg.norm(eCC),1) == 1.0),\
            "|eCC| = " + str(numpy.linalg.norm(eCC)) + "≠ 1. Aborting"
        assert (round(numpy.linalg.norm(eLL),1) == 1.0),\
            "|eLL| = " + str(numpy.linalg.norm(eLL)) + "≠ 1. Aborting"

        angle_helix = farray_angle_helix.GetTuple(k_tuple)[0]
        if (angles_in_degrees): angle_helix = angle_helix*math.pi/180
        eF = math.cos(angle_helix) * eCC + math.sin(angle_helix) * eLL
        #print "eF = " + str(eF)
        if (shuffle_vectors):
            eF *= random.choice([-1,+1])
            #print "eF = " + str(eF)
        if (use_new_definition):
            eN = eRR
            if (shuffle_vectors):
                eN *= random.choice([-1,+1])
                assert (abs(numpy.dot(eN,eRR)) > 0.999)
            eS = numpy.cross(eN, eF)
        else:
            eS = eRR
            if (shuffle_vectors): eS *= random.choice([-1,+1])
            eN = numpy.cross(eF, eS)

        assert (round(numpy.linalg.norm(eF),1) == 1.0),\
            "|eF| = " + str(numpy.linalg.norm(eF)) + "≠ 1. Aborting"
        assert (round(numpy.linalg.norm(eS),1) == 1.0),\
            "|eS| = " + str(numpy.linalg.norm(eS)) + "≠ 1. Aborting"
        assert (round(numpy.linalg.norm(eN),1) == 1.0),\
            "|eN| = " + str(numpy.linalg.norm(eN)) + "≠ 1. Aborting"

        farray_eF.InsertTuple(k_tuple, eF)
        farray_eS.InsertTuple(k_tuple, eS)
        farray_eN.InsertTuple(k_tuple, eN)

    return (farray_eF,
            farray_eS,
            farray_eN)

########################################################################

def addFiberDirections(
        ugrid,
        type_of_support="cell",
        angles_in_degrees=True,
        use_new_definition=False,
        shuffle_vectors=False,
        verbose=1):

    myVTK.myPrint(verbose, "*** addFiberDirections ***")

    if (type_of_support == "cell"):
        ugrid_data = ugrid.GetCellData()
    elif (type_of_support == "point"):
        ugrid_data = ugrid.GetPointData()

    farray_eRR = ugrid_data.GetArray("eRR")
    farray_eCC = ugrid_data.GetArray("eCC")
    farray_eLL = ugrid_data.GetArray("eLL")

    farray_angle_helix = ugrid_data.GetArray("angle_helix")

    (farray_eF,
     farray_eS,
     farray_eN) = computeFiberDirections(
        farray_eRR=farray_eRR,
        farray_eCC=farray_eCC,
        farray_eLL=farray_eLL,
        farray_angle_helix=farray_angle_helix,
        angles_in_degrees=angles_in_degrees,
        use_new_definition=use_new_definition,
        shuffle_vectors=shuffle_vectors,
        verbose=verbose-1)

    ugrid_data.AddArray(farray_eF)
    ugrid_data.AddArray(farray_eS)
    ugrid_data.AddArray(farray_eN)

    return (farray_eF,
            farray_eS,
            farray_eN)
