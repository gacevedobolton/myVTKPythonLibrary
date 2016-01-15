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

import math
import numpy

import myVTKPythonLibrary as myVTK

########################################################################

def computeHelixTransverseSheetAngles2(
        farray_eRR,
        farray_eCC,
        farray_eLL,
        farray_eF,
        farray_eS,
        farray_eN,
        use_new_definition=False,
        ref_vectors_are_material_basis=False,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeHelixTransverseSheetAngles2 ***")

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
        #print "k_tuple = " + str(k_tuple)

        eRR = numpy.array(farray_eRR.GetTuple(k_tuple))
        eCC = numpy.array(farray_eCC.GetTuple(k_tuple))
        eLL = numpy.array(farray_eLL.GetTuple(k_tuple))
        eF  = numpy.array(farray_eF.GetTuple(k_tuple))
        eS  = numpy.array(farray_eS.GetTuple(k_tuple))
        eN  = numpy.array(farray_eN.GetTuple(k_tuple))

        #print "eRR = " + str(eRR)
        #print "eCC = " + str(eCC)
        #print "eLL = " + str(eLL)
        #print "eF = " + str(eF)
        #print "eS = " + str(eS)
        #print "eN = " + str(eN)

        #print "|eRR| = " + str(numpy.linalg.norm(eRR))
        #print "|eCC| = " + str(numpy.linalg.norm(eCC))
        #print "|eLL| = " + str(numpy.linalg.norm(eLL))

        #print "eRR.eCC = " + str(numpy.dot(eRR, eCC))
        #print "eCC.eLL = " + str(numpy.dot(eCC, eLL))
        #print "eLL.eRR = " + str(numpy.dot(eLL, eRR))

        #print "|eF| = " + str(numpy.linalg.norm(eF))
        #print "|eS| = " + str(numpy.linalg.norm(eS))
        #print "|eN| = " + str(numpy.linalg.norm(eN))

        #print "eF.eS = " + str(numpy.dot(eF, eS))
        #print "eS.eN = " + str(numpy.dot(eS, eN))
        #print "eN.eF = " + str(numpy.dot(eN, eF))

        assert (round(numpy.linalg.norm(eRR),1) == 1.0),\
            "|eRR| = " + str(numpy.linalg.norm(eRR)) + " ≠ 1. Aborting."
        assert (round(numpy.linalg.norm(eCC),1) == 1.0),\
            "|eCC| = " + str(numpy.linalg.norm(eCC)) + " ≠ 1. Aborting."
        assert (round(numpy.linalg.norm(eLL),1) == 1.0),\
            "|eLL| = " + str(numpy.linalg.norm(eLL)) + " ≠ 1. Aborting."

        assert (round(numpy.dot(eRR,eCC),1) == 0.0),\
            "eRR.eCC = " + str(numpy.dot(eRR,eCC)) + " ≠ 0. Aborting."
        assert (round(numpy.dot(eCC,eLL),1) == 0.0),\
            "eCC.eLL = " + str(numpy.dot(eCC,eLL)) + " ≠ 0. Aborting."
        assert (round(numpy.dot(eLL,eRR),1) == 0.0),\
            "eLL.eRR = " + str(numpy.dot(eLL,eRR)) + " ≠ 0. Aborting."

        assert (round(numpy.linalg.det([eRR, eCC, eLL]),1) == 1.0),\
            "det([eRR, eCC, eLL]) = " + str(numpy.linalg.det([eRR, eCC, eLL])) + " ≠ 1. Aborting."

        assert (round(numpy.linalg.norm(eF),1) == 1.0),\
            "|eF| = " + str(numpy.linalg.norm(eF)) + " ≠ 1. Aborting."
        assert (round(numpy.linalg.norm(eS),1) == 1.0),\
            "|eS| = " + str(numpy.linalg.norm(eS)) + " ≠ 1. Aborting."
        assert (round(numpy.linalg.norm(eN),1) == 1.0),\
            "|eN| = " + str(numpy.linalg.norm(eN)) + " ≠ 1. Aborting."

        assert (round(numpy.dot(eF,eS),1) == 0.0),\
            "eF.eS = " + str(numpy.dot(eF,eS)) + " ≠ 0. Aborting."
        assert (round(numpy.dot(eS,eN),1) == 0.0),\
            "eS.eN = " + str(numpy.dot(eS,eN)) + " ≠ 0. Aborting."
        assert (round(numpy.dot(eN,eF),1) == 0.0),\
            "eN.eF = " + str(numpy.dot(eN,eF)) + " ≠ 0. Aborting."

        assert (round(numpy.linalg.det([eF, eS, eN]),1) == 1.0),\
            "det([eF, eS, eN]) = " + str(numpy.linalg.det([eF, eS, eN])) + " ≠ 1. Aborting."

        # reference basis
        if (ref_vectors_are_material_basis):
            ref = numpy.array([+eRR, +eCC, +eLL])
        else:
            if (use_new_definition):
                ref = numpy.array([+eCC, +eLL, +eRR])
            else:
                ref = numpy.array([+eCC, +eRR, -eLL])

        # material basis
        if (abs(numpy.dot(eF, ref[0])) > 1e-3):
            eF = math.copysign(1., numpy.dot(eF, ref[0])) * eF
        if (abs(numpy.dot(eS, ref[1])) > 1e-3):
            eS = math.copysign(1., numpy.dot(eS, ref[1])) * eS
        eN = numpy.cross(eF, eS)
        base = numpy.array([eF, eS, eN])

        assert (round(numpy.linalg.det(base),1) == 1.0),\
            "det(base) = " + str(numpy.linalg.det(base)) + " ≠ 1. Aborting."

        #print "b0.r0 = " + str(numpy.dot(base[0], ref[0]))
        #print "b0.r1 = " + str(numpy.dot(base[0], ref[1]))
        #print "b0.r2 = " + str(numpy.dot(base[0], ref[2]))
        #print "b1.r0 = " + str(numpy.dot(base[1], ref[0]))
        #print "b1.r1 = " + str(numpy.dot(base[1], ref[1]))
        #print "b1.r2 = " + str(numpy.dot(base[1], ref[2]))
        #print "b2.r0 = " + str(numpy.dot(base[2], ref[0]))
        #print "b2.r1 = " + str(numpy.dot(base[2], ref[1]))
        #print "b2.r2 = " + str(numpy.dot(base[2], ref[2]))

        if (use_new_definition):
            sheet = math.atan2(numpy.dot(base[1], ref[2]), numpy.dot(base[2], ref[2]))
            #print "sheet = " + str(sheet) + " = " + str(sheet*180/math.pi)
            #sheet = (sheet+math.pi/2)%math.pi - math.pi/2
            #print "sheet = " + str(sheet) + " = " + str(sheet*180/math.pi)
            C = math.cos(sheet)
            S = math.sin(sheet)
            R_sheet = numpy.array([[1,  0, 0],\
                                   [0,  C, S],\
                                   [0, -S, C]])
            base = numpy.dot(numpy.transpose(R_sheet), base)
            #print base
            if (numpy.dot(base[1], ref[1]) < 0):
                R_sheet = numpy.array([[1,  0,  0],\
                                       [0, -1,  0],\
                                       [0,  0, -1]])
                base = numpy.dot(numpy.transpose(R_sheet), base)
        else:
            sheet = math.atan2(-numpy.dot(base[2], ref[1]), numpy.dot(base[1], ref[1]))
            #print "sheet = " + str(sheet) + " = " + str(sheet*180/math.pi)
            #sheet = (sheet+math.pi/2)%math.pi - math.pi/2
            #print "sheet = " + str(sheet) + " = " + str(sheet*180/math.pi)
            C = math.cos(sheet)
            S = math.sin(sheet)
            R_sheet = numpy.array([[1,  0, 0],\
                                   [0,  C, S],\
                                   [0, -S, C]])
            base = numpy.dot(numpy.transpose(R_sheet), base)
            #print base
            if (numpy.dot(base[2], ref[2]) < 0):
                R_sheet = numpy.array([[1,  0,  0],\
                                       [0, -1,  0],\
                                       [0,  0, -1]])
                base = numpy.dot(numpy.transpose(R_sheet), base)

        #print "b0.r0 = " + str(numpy.dot(base[0], ref[0]))
        #print "b0.r1 = " + str(numpy.dot(base[0], ref[1]))
        #print "b0.r2 = " + str(numpy.dot(base[0], ref[2]))
        #print "b1.r0 = " + str(numpy.dot(base[1], ref[0]))
        #print "b1.r1 = " + str(numpy.dot(base[1], ref[1]))
        #print "b1.r2 = " + str(numpy.dot(base[1], ref[2]))
        #print "b2.r0 = " + str(numpy.dot(base[2], ref[0]))
        #print "b2.r1 = " + str(numpy.dot(base[2], ref[1]))
        #print "b2.r2 = " + str(numpy.dot(base[2], ref[2]))

        if (use_new_definition):
            trans = math.atan2(numpy.dot(base[2], ref[0]), numpy.dot(base[0], ref[0]))
            #print "trans = " + str(trans) + " = " + str(trans*180/math.pi)
            #assert (math.atan2(numpy.dot(base[2], ref[1]), numpy.dot(base[0], ref[1])) == trans),\
                #"atan2(b2 . r1, b0 . r1) = " + str(math.atan2(numpy.dot(base[2], ref[1]), numpy.dot(base[0], ref[1]))) + " ≠ trans. Aborting."
            #assert (math.atan2(-numpy.dot(base[0], ref[2]), numpy.dot(base[2], ref[2])) == trans),\
                #"atan2(-b0 . r2, b2 . r2) = " + str(math.atan2(-numpy.dot(base[0], ref[2]), numpy.dot(base[2], ref[2]))) + " ≠ trans. Aborting."
            #trans = (trans+math.pi/2)%math.pi - math.pi/2
            #print "trans = " + str(trans) + " = " + str(trans*180/math.pi)
            C = math.cos(trans)
            S = math.sin(trans)
            R_trans = numpy.array([[ C, 0,-S],\
                                   [ 0, 1, 0],\
                                   [ S, 0, C]])
            base = numpy.dot(numpy.transpose(R_trans), base)
            #print base
            if (numpy.dot(base[2], ref[2]) < 0):
                R_trans = numpy.array([[-1,  0,  0],\
                                       [ 0,  1,  0],\
                                       [ 0,  0, -1]])
                base = numpy.dot(numpy.transpose(R_trans), base)
        else:
            trans = math.atan2(-numpy.dot(base[1], ref[0]), numpy.dot(base[0], ref[0]))
            #print "trans = " + str(trans) + " = " + str(trans*180/math.pi)
            #assert (math.atan2(-numpy.dot(base[1], ref[2]), numpy.dot(base[0], ref[2])) == trans),\
                #"atan2(b1 . r2, b0 . r2) = " + str(math.atan2(numpy.dot(base[1], ref[2]), numpy.dot(base[0], ref[2]))) + " ≠ trans. Aborting."
            #assert (-math.atan2(-numpy.dot(base[0], ref[1]), numpy.dot(base[1], ref[1])) == trans),\
                #"atan2(-b0 . r1, b1 . r1) = " + str(math.atan2(-numpy.dot(base[0], ref[1]), numpy.dot(base[1], ref[1]))) + " ≠ trans. Aborting."
            #trans = (trans+math.pi/2)%math.pi - math.pi/2
            #print "trans = " + str(trans) + " = " + str(trans*180/math.pi)
            C = math.cos(trans)
            S = math.sin(trans)
            R_trans = numpy.array([[ C, S, 0],\
                                   [-S, C, 0],\
                                   [ 0, 0, 1]])
            base = numpy.dot(numpy.transpose(R_trans), base)
            #print base
            if (numpy.dot(base[1], ref[1]) < 0):
                R_trans = numpy.array([[-1,  0,  0],\
                                       [ 0, -1,  0],\
                                       [ 0,  0,  1]])
                base = numpy.dot(numpy.transpose(R_trans), base)

        #print "b0.r0 = " + str(numpy.dot(base[0], ref[0]))
        #print "b0.r1 = " + str(numpy.dot(base[0], ref[1]))
        #print "b0.r2 = " + str(numpy.dot(base[0], ref[2]))
        #print "b1.r0 = " + str(numpy.dot(base[1], ref[0]))
        #print "b1.r1 = " + str(numpy.dot(base[1], ref[1]))
        #print "b1.r2 = " + str(numpy.dot(base[1], ref[2]))
        #print "b2.r0 = " + str(numpy.dot(base[2], ref[0]))
        #print "b2.r1 = " + str(numpy.dot(base[2], ref[1]))
        #print "b2.r2 = " + str(numpy.dot(base[2], ref[2]))

        if (use_new_definition):
            helix = math.atan2(numpy.dot(base[0], ref[1]), numpy.dot(base[1], ref[1]))
            #print "helix = " + str(helix) + " = " + str(helix*180/math.pi)
            #assert (numpy.isclose(math.atan2(-numpy.dot(base[1], ref[0]), numpy.dot(base[0], ref[0])), helix, atol=1e-3)),\
                #"atan2(-b1 . r0, b0 . r0) = " + str(math.atan2(-numpy.dot(base[1], ref[0]), numpy.dot(base[0], ref[0]))) + " ≠ helix. Aborting."
            #helix = (helix+math.pi/2)%math.pi - math.pi/2
            #print "helix = " + str(helix) + " = " + str(helix*180/math.pi)
            C = math.cos(helix)
            S = math.sin(helix)
            R_helix = numpy.array([[ C, S, 0],\
                                   [-S, C, 0],\
                                   [ 0, 0, 1]])
            base = numpy.dot(numpy.transpose(R_helix), base)
            #print base
            if (numpy.dot(base[0], ref[0]) < 0):
                R_helix = numpy.array([[-1,  0,  0],\
                                       [ 0, -1,  0],\
                                       [ 0,  0,  1]])
                base = numpy.dot(numpy.transpose(R_helix), base)
        else:
            helix = math.atan2(-numpy.dot(base[0], ref[2]), numpy.dot(base[2], ref[2]))
            #print "helix = " + str(helix) + " = " + str(helix*180/math.pi)
            #assert (numpy.isclose(math.atan2(numpy.dot(base[2], ref[0]), numpy.dot(base[0], ref[0])), helix, atol=1e-3)),\
                #"atan2(-b2 . r0, b0 . r0) = " + str(math.atan2(-numpy.dot(base[2], ref[0]), numpy.dot(base[0], ref[0]))) + " ≠ helix. Aborting."
            #helix = (helix+math.pi/2)%math.pi - math.pi/2
            #print "helix = " + str(helix) + " = " + str(helix*180/math.pi)
            C = math.cos(helix)
            S = math.sin(helix)
            R_helix = numpy.array([[ C, 0,-S],\
                                   [ 0, 1, 0],\
                                   [ S, 0, C]])
            base = numpy.dot(numpy.transpose(R_helix), base)
            #print base
            if (numpy.dot(base[0], ref[0]) < 0):
                R_helix = numpy.array([[-1,  0,  0],\
                                       [ 0,  1,  0],\
                                       [ 0,  0, -1]])
                base = numpy.dot(numpy.transpose(R_helix), base)

        #print "b0.r0 = " + str(numpy.dot(base[0], ref[0]))
        #print "b0.r1 = " + str(numpy.dot(base[0], ref[1]))
        #print "b0.r2 = " + str(numpy.dot(base[0], ref[2]))
        #print "b1.r0 = " + str(numpy.dot(base[1], ref[0]))
        #print "b1.r1 = " + str(numpy.dot(base[1], ref[1]))
        #print "b1.r2 = " + str(numpy.dot(base[1], ref[2]))
        #print "b2.r0 = " + str(numpy.dot(base[2], ref[0]))
        #print "b2.r1 = " + str(numpy.dot(base[2], ref[1]))
        #print "b2.r2 = " + str(numpy.dot(base[2], ref[2]))

        assert (round(numpy.dot(base[0],ref[0]),1) == 1.0),\
            "b0.r0 = " + str(numpy.dot(base[0],ref[0])) + " ≠ 1. Aborting."
        assert (round(numpy.dot(base[1],ref[1]),1) == 1.0),\
            "b1.r1 = " + str(numpy.dot(base[1],ref[1])) + " ≠ 1. Aborting."
        assert (round(numpy.dot(base[2],ref[2]),1) == 1.0),\
            "b2.r2 = " + str(numpy.dot(base[2],ref[2])) + " ≠ 1. Aborting."

        helix = (helix+math.pi/2)%math.pi - math.pi/2
        trans = (trans+math.pi/2)%math.pi - math.pi/2
        sheet = (sheet+math.pi/2)%math.pi - math.pi/2

        helix *= 180/math.pi
        trans *= 180/math.pi
        sheet *= 180/math.pi

        farray_angle_helix.InsertTuple(k_tuple, [helix])
        farray_angle_trans.InsertTuple(k_tuple, [trans])
        farray_angle_sheet.InsertTuple(k_tuple, [sheet])

    return (farray_angle_helix,
            farray_angle_trans,
            farray_angle_sheet)

########################################################################

def addHelixTransverseSheetAngles2(
        ugrid,
        type_of_support="cell",
        use_new_definition=False,
        ref_vectors_are_material_basis=False,
        verbose=1):

    myVTK.myPrint(verbose, "*** addHelixTransverseSheetAngles2 ***")

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
     farray_angle_sheet) = computeHelixTransverseSheetAngles2(
        farray_eRR=farray_eRR,
        farray_eCC=farray_eCC,
        farray_eLL=farray_eLL,
        farray_eF=farray_eF,
        farray_eS=farray_eS,
        farray_eN=farray_eN,
        use_new_definition=use_new_definition,
        ref_vectors_are_material_basis=ref_vectors_are_material_basis,
        verbose=verbose-1)

    ugrid_data.AddArray(farray_angle_helix)
    ugrid_data.AddArray(farray_angle_trans)
    ugrid_data.AddArray(farray_angle_sheet)

    return (farray_angle_helix,
            farray_angle_trans,
            farray_angle_sheet)
