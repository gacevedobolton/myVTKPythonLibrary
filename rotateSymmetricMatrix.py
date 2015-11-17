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
import vtk

import myVTKPythonLibrary as myVTK
from mat_vec_tools import *

########################################################################

def rotateSymmetricMatrix(
        old_array,
        old_array_storage="vec",
        in_vecs=None,
        out_vecs=None,
        verbose=1):

    myVTK.myPrint(verbose, "*** rotateSymmetricMatrix ***")

    n_tuples = old_array.GetNumberOfTuples()
    new_array = myVTK.createFloatArray("", 6, n_tuples)

    for k_tuple in xrange(n_tuples):
        old_matrix = old_array.GetTuple(k_tuple)
        if (old_array_storage == "vec"):
            old_matrix = vec_col_to_mat_sym(old_matrix)
        elif (old_array_storage == "Cmat"):
            old_matrix = numpy.reshape(old_matrix, (3,3), order='C')
        elif (old_array_storage == "Fmat"):
            old_matrix = numpy.reshape(old_matrix, (3,3), order='F')

        if (in_vecs == None):
            in_R = numpy.eye(3)
        else:
            in_R = numpy.transpose(numpy.array([in_vecs[0].GetTuple(k_tuple),
                                                in_vecs[1].GetTuple(k_tuple),
                                                in_vecs[2].GetTuple(k_tuple)]))

        if (out_vecs == None):
            out_R = numpy.eye(3)
        else:
            out_R = numpy.transpose(numpy.array([out_vecs[0].GetTuple(k_tuple),
                                                 out_vecs[1].GetTuple(k_tuple),
                                                 out_vecs[2].GetTuple(k_tuple)]))

        R = numpy.dot(numpy.transpose(in_R), out_R)

        new_matrix = numpy.dot(numpy.dot(numpy.transpose(R), old_matrix), R)

        if (old_array_storage == "vec"):
            new_matrix = mat_sym_to_vec_col(new_matrix)
        elif (old_array_storage == "Cmat"):
            new_matrix = numpy.reshape(new_matrix, 9, order='C')
        elif (old_array_storage == "Fmat"):
            new_matrix = numpy.reshape(new_matrix, 9, order='F')

        new_array.InsertTuple(k_tuple, new_matrix)

    return new_array
