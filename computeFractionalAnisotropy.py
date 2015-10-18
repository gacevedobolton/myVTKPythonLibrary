#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import myVTKPythonLibrary as myVTK

########################################################################

def computeFractionalAnisotropy(
        farray_e1,
        farray_e2,
        farray_e3,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeFractionalAnisotropy ***")

    n_tuples = farray_e1.GetNumberOfTuples()

    farray_FA   = myVTK.createFloatArray("FA"   , 1, n_tuples)
    farray_FA12 = myVTK.createFloatArray("FA_12", 1, n_tuples)
    farray_FA23 = myVTK.createFloatArray("FA_23", 1, n_tuples)

    for k_tuple in xrange(n_tuples):
        e1 = farray_e1.GetTuple(k_tuple)[0]
        e2 = farray_e2.GetTuple(k_tuple)[0]
        e3 = farray_e3.GetTuple(k_tuple)[0]
        FA   = ((e1-e2)**2+(e1-e3)**2+(e2-e3)**2)**(0.5) / (2*(e1**2+e2**2+e3**2))**(0.5)
        FA12 = ((e1-e2)**2)**(0.5) / (e1**2+e2**2)**(0.5)
        FA23 = ((e2-e3)**2)**(0.5) / (e2**2+e3**2)**(0.5)

        farray_FA.InsertTuple(k_tuple, [FA])
        farray_FA12.InsertTuple(k_tuple, [FA12])
        farray_FA23.InsertTuple(k_tuple, [FA23])

    return (farray_FA,
            farray_FA12,
            farray_FA23)

########################################################################

def addFractionalAnisotropy(
        ugrid,
        field_name,
        type_of_support="cell",
        verbose=1):

    myVTK.myPrint(verbose, "*** addFractionalAnisotropy ***")

    if (type_of_support == "cell"):
        data = ugrid.GetCellData()
    elif (type_of_support == "point"):
        data = ugrid.GetPointData()

    farray_e1 = data.GetArray(field_name+"_Lmax")
    farray_e2 = data.GetArray(field_name+"_Lmid")
    farray_e3 = data.GetArray(field_name+"_Lmin")

    (farray_FA,
    farray_FA12,
    farray_FA23) = computeFractionalAnisotropy(
        farray_e1=farray_e1,
        farray_e2=farray_e2,
        farray_e3=farray_e3,
        verbose=verbose-1)

    farray_FA.SetName(field_name+"_FA")
    farray_FA12.SetName(field_name+"_FA_12")
    farray_FA23.SetName(field_name+"_FA_23")

    data.AddArray(farray_FA)
    data.AddArray(farray_FA12)
    data.AddArray(farray_FA23)

    return (farray_FA,
            farray_FA12,
            farray_FA23)
