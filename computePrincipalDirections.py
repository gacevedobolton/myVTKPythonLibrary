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

import math
import numpy

import myVTKPythonLibrary as myVTK
from mat_vec_tools import *

########################################################################

def computePrincipalDirections(
        field,
        field_storage="vec",
        orient=0,
        farray_eRR=None,
        farray_eCC=None,
        farray_eLL=None,
        verbose=1):

    myVTK.myPrint(verbose, "*** computePrincipalDirections ***")

    assert (field_storage in ["vec", "Cmat", "Fmat"]), "\"field_storage\" must be \"vec\", \"Cmat\" or \"Fmat\". Aborting."

    n_tuples = field.GetNumberOfTuples()

    farray_Lmin = myVTK.createFloatArray('Lmin', 1, n_tuples)
    farray_Lmid = myVTK.createFloatArray('Lmid', 1, n_tuples)
    farray_Lmax = myVTK.createFloatArray('Lmax', 1, n_tuples)

    farray_Vmin = myVTK.createFloatArray('Vmin', 3, n_tuples)
    farray_Vmid = myVTK.createFloatArray('Vmid', 3, n_tuples)
    farray_Vmax = myVTK.createFloatArray('Vmax', 3, n_tuples)

    for k_tuple in xrange(n_tuples):
        #print "k_tuple: " + str(k_tuple)
        matrix = field.GetTuple(k_tuple)
        if (field_storage == "vec"):
            matrix = vec_col_to_mat_sym(matrix)
        elif (field_storage == "Cmat"):
            matrix = numpy.reshape(matrix, (3,3), order='C')
        elif (field_storage == "Fmat"):
            matrix = numpy.reshape(matrix, (3,3), order='F')

        if (numpy.linalg.norm(matrix) > 1e-6):
            #if (verbose): print indent(verbose) + 'k_tuple =', k_tuple

            val, vec = numpy.linalg.eig(matrix)
            #if (verbose): print indent(verbose) + 'val =', val
            #if (verbose): print indent(verbose) + 'vec =', vec
            #if (verbose): print indent(verbose) + 'det =', numpy.linalg.det(vec)
            idx = val.argsort()
            val = val[idx]
            vec = vec[:,idx]
            #if (verbose): print indent(verbose) + 'val =', val
            #if (verbose): print indent(verbose) + 'vec =', vec
            #if (verbose): print indent(verbose) + 'det =', numpy.linalg.det(vec)

            matrix_Lmin = [val[0]]
            matrix_Lmid = [val[1]]
            matrix_Lmax = [val[2]]

            matrix_Vmax = vec[:,2]
            matrix_Vmid = vec[:,1]
            if (orient):
                eCC = numpy.array(farray_eCC.GetTuple(k_tuple))
                eLL = numpy.array(farray_eLL.GetTuple(k_tuple))
                matrix_Vmax = math.copysign(1, numpy.dot(matrix_Vmax, eCC)) * matrix_Vmax
                matrix_Vmid = math.copysign(1, numpy.dot(matrix_Vmid, eLL)) * matrix_Vmid
            matrix_Vmin = numpy.cross(matrix_Vmax, matrix_Vmid)
        else:
            matrix_Lmin = [0.]
            matrix_Lmid = [0.]
            matrix_Lmax = [0.]
            matrix_Vmin = [0.]*3
            matrix_Vmid = [0.]*3
            matrix_Vmax = [0.]*3

        farray_Lmin.InsertTuple(k_tuple, matrix_Lmin)
        farray_Lmid.InsertTuple(k_tuple, matrix_Lmid)
        farray_Lmax.InsertTuple(k_tuple, matrix_Lmax)
        farray_Vmin.InsertTuple(k_tuple, matrix_Vmin)
        farray_Vmid.InsertTuple(k_tuple, matrix_Vmid)
        farray_Vmax.InsertTuple(k_tuple, matrix_Vmax)

    return (farray_Lmin,
            farray_Lmid,
            farray_Lmax,
            farray_Vmin,
            farray_Vmid,
            farray_Vmax)

########################################################################

def addPrincipalDirections(
        ugrid,
        field_name,
        field_support="cell",
        field_storage="vec",
        orient=0,
        verbose=1):

    myVTK.myPrint(verbose, "*** addPrincipalDirections ***")

    assert (field_support in ["point", "cell"]), "\"field_support\" must be \"point\" or \"cell\". Aborting."
    assert (field_storage in ["vec", "Cmat", "Fmat"]), "\"field_storage\" must be \"vec\", \"Cmat\" or \"Fmat\". Aborting."

    if   (field_support == "cell" ): ugrid_data = ugrid.GetCellData()
    elif (field_support == "point"): ugrid_data = ugrid.GetPointData()

    field = ugrid_data.GetArray(field_name)

    if (orient):
        (farray_Lmin,
        farray_Lmid,
        farray_Lmax,
        farray_Vmin,
        farray_Vmid,
        farray_Vmax) = computePrincipalDirections(
            field=field,
            field_storage=field_storage,
            orient=orient,
            farray_eRR=ugrid_data.GetArray("eRR"),
            farray_eCC=ugrid_data.GetArray("eCC"),
            farray_eLL=ugrid_data.GetArray("eLL"),
            verbose=verbose-1)
    else:
        (farray_Lmin,
        farray_Lmid,
        farray_Lmax,
        farray_Vmin,
        farray_Vmid,
        farray_Vmax) = computePrincipalDirections(
            field=field,
            field_storage=field_storage,
            orient=orient,
            verbose=verbose-1)

    farray_Lmin.SetName(field_name+"_Lmin")
    farray_Lmid.SetName(field_name+"_Lmid")
    farray_Lmax.SetName(field_name+"_Lmax")
    farray_Vmin.SetName(field_name+"_Vmin")
    farray_Vmid.SetName(field_name+"_Vmid")
    farray_Vmax.SetName(field_name+"_Vmax")

    ugrid_data.AddArray(farray_Lmin)
    ugrid_data.AddArray(farray_Lmid)
    ugrid_data.AddArray(farray_Lmax)
    ugrid_data.AddArray(farray_Vmin)
    ugrid_data.AddArray(farray_Vmid)
    ugrid_data.AddArray(farray_Vmax)

    return (farray_Lmin,
            farray_Lmid,
            farray_Lmax,
            farray_Vmin,
            farray_Vmid,
            farray_Vmax)
