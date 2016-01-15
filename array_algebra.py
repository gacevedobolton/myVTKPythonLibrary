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

import sys
import math
import random
import numpy

import myVTKPythonLibrary as myVTK

########################################################################

def addArrays(
    array1,
    array2,
    array3=None,
    verbose=0):

    myVTK.myPrint(verbose, "*** addArrays ***")

    n_components = array1.GetNumberOfComponents()
    assert (array2.GetNumberOfComponents() == n_components)
    if (verbose >= 2): print "n_components = " + str(n_components)

    n_tuples = array1.GetNumberOfTuples()
    assert (array2.GetNumberOfTuples() == n_tuples)
    if (verbose >= 2): print "n_tuples = " + str(n_tuples)

    array_type = type(array1.GetTuple(0)[0])
    assert (array_type in [int, float])
    assert (type(array2.GetTuple(0)[0]) is array_type)
    if (verbose >= 2): print "array_type = " + str(array_type)

    if (array3 == None):
        array3 = myVTK.createArray(
            array_name="",
            n_components=n_components,
            n_tuples=n_tuples,
            array_type=array_type)
    else:
        assert (array3.GetNumberOfComponents() == n_components)
        assert (array3.GetNumberOfTuples() == n_tuples)
        assert (type(array3.GetTuple(0)[0]) is array_type)

    for k_tuple in xrange(n_tuples):
        if (verbose >= 2): print "k_tuple = " + str(k_tuple)

        array3.SetTuple(
            k_tuple,
            numpy.array(array1.GetTuple(k_tuple)) + numpy.array(array2.GetTuple(k_tuple)))

    return array3

########################################################################

def subArrays(
    array1,
    array2,
    array3=None,
    verbose=0):

    myVTK.myPrint(verbose, "*** subArrays ***")

    n_components = array1.GetNumberOfComponents()
    assert (array2.GetNumberOfComponents() == n_components)

    n_tuples = array1.GetNumberOfTuples()
    assert (array2.GetNumberOfTuples() == n_tuples)

    array_type = type(array1.GetTuple(0)[0])
    assert (array_type in [int, float])
    assert (type(array2.GetTuple(0)[0]) is array_type)

    if (array3 == None):
        array3 = myVTK.createArray(
            array_name="",
            n_components=n_components,
            n_tuples=n_tuples,
            array_type=array_type)
    else:
        assert (array3.GetNumberOfComponents() == n_components)
        assert (array3.GetNumberOfTuples() == n_tuples)
        assert (type(array3.GetTuple(0)[0]) is array_type)

    for k_tuple in xrange(n_tuples):
        array3.SetTuple(
            k_tuple,
            numpy.array(array1.GetTuple(k_tuple)) - numpy.array(array2.GetTuple(k_tuple)))

    return array3

########################################################################

def mulArrays(
    array1,
    array2,
    array3=None,
    verbose=0):

    myVTK.myPrint(verbose, "*** mulArrays ***")

    n_components = array1.GetNumberOfComponents()
    assert (array2.GetNumberOfComponents() == n_components)

    n_tuples = array1.GetNumberOfTuples()
    assert (array2.GetNumberOfTuples() == n_tuples)

    array_type = type(array1.GetTuple(0)[0])
    assert (array_type in [int, float])
    assert (type(array2.GetTuple(0)[0]) is array_type)

    if (array3 == None):
        array3 = myVTK.createArray(
            array_name="",
            n_components=n_components,
            n_tuples=n_tuples,
            array_type=array_type)
    else:
        assert (array3.GetNumberOfComponents() == n_components)
        assert (array3.GetNumberOfTuples() == n_tuples)
        assert (type(array3.GetTuple(0)[0]) is array_type)

    for k_tuple in xrange(n_tuples):
        array3.SetTuple(
            k_tuple,
            numpy.array(array1.GetTuple(k_tuple)) * numpy.array(array2.GetTuple(k_tuple)))

    return array3
