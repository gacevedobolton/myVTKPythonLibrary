#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import numpy

import myVTKPythonLibrary as myVTK

########################################################################

def readDynaDeformationGradients(
        mesh,
        hystory_files_basename,
        array_name,
        verbose=1):

    myVTK.myPrint(verbose, "*** readDynaDeformationGradients ***")

    n_cells = mesh.GetNumberOfCells()

    history_files_names = [hystory_files_basename + '.history#' + str(num) for num in xrange(11,20)]

    F_list = [[0. for k_component in xrange(9)] for k_cell in xrange(n_cells)]

    for k_component in xrange(9):
        history_file = open(history_files_names[k_component], 'r')
        for line in history_file:
            if line.startswith('*') or line.startswith('$'): continue
            line = line.split()
            F_list[int(line[0])-1][k_component] = float(line[1])
        history_file.close()

    F_array = myVTK.createFloatArray(array_name, 9, n_cells)

    for k_cell in xrange(n_cells):
        F_array.InsertTuple(k_cell, F_list[k_cell])

    myVTK.myPrint(verbose, "n_tuples = " + str(F_array.GetNumberOfTuples()))

    mesh.GetCellData().AddArray(F_array)
