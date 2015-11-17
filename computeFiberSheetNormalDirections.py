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

import sys
import math
import random
import numpy

import myVTKPythonLibrary as myVTK

########################################################################

def computeFiberSheetNormalDirections(
        farray_eRR,
        farray_eCC,
        farray_eLL,
        farray_helix,
        farray_trans,
        farray_sheet,
        angles_in_degrees=True,
        use_new_definition=False,
        shuffle_vectors=False,
        verbose=True):

    myVTK.myPrint(verbose, "*** computeFiberSheetNormalDirections ***")

    n_tuples = farray_helix.GetNumberOfTuples()

    farray_eF = myVTK.createFloatArray("eF", 3, n_tuples)
    farray_eS = myVTK.createFloatArray("eS", 3, n_tuples)
    farray_eN = myVTK.createFloatArray("eN", 3, n_tuples)

    for k_tuple in xrange(n_tuples):
        eRR = numpy.array(farray_eRR.GetTuple(k_tuple))
        eCC = numpy.array(farray_eCC.GetTuple(k_tuple))
        eLL = numpy.array(farray_eLL.GetTuple(k_tuple))

        if (use_new_definition):
            base = numpy.array([eCC,\
                                eLL,\
                                eRR])

            helix = farray_helix.GetTuple(k_tuple)[0]
            if (angles_in_degrees): helix *= math.pi/180
            C = math.cos(helix)
            S = math.sin(helix)
            R_helix = numpy.array([[ C, S, 0],\
                                   [-S, C, 0],\
                                   [ 0, 0, 1]])
            base = numpy.dot(R_helix, base)

            trans = farray_trans.GetTuple(k_tuple)[0]
            if (angles_in_degrees): trans *= math.pi/180
            C = math.cos(trans)
            S = math.sin(trans)
            R_trans = numpy.array([[ C, 0,-S],\
                                   [ 0, 1, 0],\
                                   [ S, 0, C]])
            base = numpy.dot(R_trans, base)

            sheet = farray_sheet.GetTuple(k_tuple)[0]
            if (angles_in_degrees): sheet *= math.pi/180
            C = math.cos(sheet)
            S = math.sin(sheet)
            R_sheet = numpy.array([[1,  0, 0],\
                                   [0,  C, S],\
                                   [0, -S, C]])
            base = numpy.dot(R_sheet, base)
        else:
            base = numpy.array([+eCC,\
                                +eRR,\
                                -eLL])

            helix = farray_helix.GetTuple(k_tuple)[0]
            if (angles_in_degrees): helix *= math.pi/180
            C = math.cos(helix)
            S = math.sin(helix)
            R_helix = numpy.array([[ C, 0,-S],\
                                   [ 0, 1, 0],\
                                   [ S, 0, C]])
            base = numpy.dot(R_helix, base)

            trans = farray_trans.GetTuple(k_tuple)[0]
            if (angles_in_degrees): trans *= math.pi/180
            C = math.cos(trans)
            S = math.sin(trans)
            R_trans = numpy.array([[ C, S, 0],\
                                   [-S, C, 0],\
                                   [ 0, 0, 1]])
            base = numpy.dot(R_trans, base)

            sheet = farray_sheet.GetTuple(k_tuple)[0]
            if (angles_in_degrees): sheet *= math.pi/180
            C = math.cos(sheet)
            S = math.sin(sheet)
            R_sheet = numpy.array([[1,  0, 0],\
                                   [0,  C, S],\
                                   [0, -S, C]])
            base = numpy.dot(R_sheet, base)

        eF = base[0]
        eS = base[1]
        eN = base[2]

        if (shuffle_vectors):
            eF *= random.choice([-1,+1])
            eS *= random.choice([-1,+1])
            eN = numpy.cross(eF, eS)

        farray_eF.InsertTuple(k_tuple, eF)
        farray_eS.InsertTuple(k_tuple, eS)
        farray_eN.InsertTuple(k_tuple, eN)

    return (farray_eF,
            farray_eS,
            farray_eN)

########################################################################

def addFiberSheetNormalDirections(
        ugrid,
        type_of_support="cell",
        angles_in_degrees=True,
        use_new_definition=False,
        shuffle_vectors=False,
        verbose=1):

    myVTK.myPrint(verbose, "*** addFiberSheetNormalDirections ***")

    if (type_of_support == "cell"):
        ugrid_data = ugrid.GetCellData()
    elif (type_of_support == "point"):
        ugrid_data = ugrid.GetPointData()

    farray_eRR = ugrid_data.GetArray("eRR")
    farray_eCC = ugrid_data.GetArray("eCC")
    farray_eLL = ugrid_data.GetArray("eLL")

    farray_helix = ugrid_data.GetArray("angle_helix")
    farray_trans = ugrid_data.GetArray("angle_trans")
    farray_sheet = ugrid_data.GetArray("angle_sheet")

    (farray_eF,
     farray_eS,
     farray_eN) = computeFiberSheetNormalDirections(
        farray_eRR=farray_eRR,
        farray_eCC=farray_eCC,
        farray_eLL=farray_eLL,
        farray_helix=farray_helix,
        farray_trans=farray_trans,
        farray_sheet=farray_sheet,
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
