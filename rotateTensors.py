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

import numpy

import myVTKPythonLibrary as myVTK

########################################################################

def rotateTensors(
        old_tensors,
        R,
        storage="vec",
        verbose=1):

    myVTK.myPrint(verbose, "*** rotateTensors ***")

    assert (storage in ["vec", "Cmat", "Fmat"]), "\"storage\" must be \"vec\", \"Cmat\" or \"Fmat\". Aborting."

    n_tuples = old_tensors.GetNumberOfTuples()
    n_components = old_tensors.GetNumberOfComponents()

    new_tensors = myVTK.createFloatArray("tensors", n_components, n_tuples)

    for k_tuple in xrange(n_tuples):
        old_tensor = old_tensors.GetTuple(k_tuple)
        if (storage == "vec"):
            old_tensor = vec_col_to_mat_sym(old_tensor)
        elif (storage == "Cmat"):
            old_tensor = numpy.reshape(old_tensor, (3,3), order='C')
        elif (storage == "Fmat"):
            old_tensor = numpy.reshape(old_tensor, (3,3), order='F')

        new_tensor = numpy.dot(R, numpy.dot(old_tensor, numpy.transpose(R)))

        if (storage == "vec"):
            new_tensor = mat_sym_to_vec_col(new_tensor)
        elif (storage == "Cmat"):
            new_tensor = numpy.reshape(new_tensor, 9, order='C')
        elif (storage == "Fmat"):
            new_tensor = numpy.reshape(new_tensor, 9, order='F')

        new_tensors.InsertTuple(k_tuple, new_tensor)

    return new_tensors
