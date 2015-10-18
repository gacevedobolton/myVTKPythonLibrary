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

def readAbaqusDeformationGradientsFromDAT(
        data_filename,
        verbose=1):

    myVTK.myPrint(verbose, "*** readAbaqusDeformationGradientsFromDAT: " + data_filename + " ***")

    farray_F = myVTK.createFloatArray("F", 9)

    data_file = open(data_filename, 'r')
    context = ""
    k_cell = 0
    for line in data_file:
        if (context == "reading deformation gradients"):
            #print line
            if ("MAXIMUM" in line):
                context = ""
                continue
            if ("OR" in line):
                splitted_line = line.split()
                assert (int(splitted_line[0]) == k_cell+1), "Wrong element number. Aborting."
                F_list = [float(splitted_line[ 3]), float(splitted_line[ 6]), float(splitted_line[7]),
                          float(splitted_line[ 9]), float(splitted_line[ 4]), float(splitted_line[8]),
                          float(splitted_line[10]), float(splitted_line[11]), float(splitted_line[5])]
                farray_F.InsertNextTuple(F_list)
                k_cell += 1

        if (line == "    ELEMENT  PT FOOT-       DG11        DG22        DG33        DG12        DG13        DG23        DG21        DG31        DG32    \n"):
            context = "reading deformation gradients"

    data_file.close()

    myVTK.myPrint(verbose, "n_tuples = " + str(farray_F.GetNumberOfTuples()))

    return farray_F





