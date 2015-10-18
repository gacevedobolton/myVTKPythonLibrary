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

def createArray(
        array_name,
        n_components=1,
        n_tuples=0,
        array_type="float",
        verbose=1):

    assert (type(array_type) in [type, str]), "array_type must be a type or a str. Aborting."

    if (type(array_type) is type):
        assert (array_type in [float, int]), "if a type, array_type must be equal to float or int. Aborting."

        if (array_type == float):
            return myVTK.createFloatArray(
                       name=array_name,
                       n_components=n_components,
                       n_tuples=n_tuples,
                       verbose=verbose-1)
        elif (array_type == int):
            return myVTK.createIntArray(
                       name=array_name,
                       n_components=n_components,
                       n_tuples=n_tuples,
                       verbose=verbose-1)
    elif (type(array_type) is str):
        assert (array_type in ["double", "float", "int", "short"]), "if a str, array_type must be equal to double, float, int or short. Aborting."

        if (array_type == "float") or (array_type == "double"):
            return myVTK.createFloatArray(
                       name=array_name,
                       n_components=n_components,
                       n_tuples=n_tuples,
                       verbose=verbose-1)
        elif (array_type == "int"):
            return myVTK.createIntArray(
                       name=array_name,
                       n_components=n_components,
                       n_tuples=n_tuples,
                       verbose=verbose-1)
        elif (array_type == "short"):
            return myVTK.createShortArray(
                       name=array_name,
                       n_components=n_components,
                       n_tuples=n_tuples,
                       verbose=verbose-1)
