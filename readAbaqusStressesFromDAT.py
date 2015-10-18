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

def readAbaqusStressFromDAT(
        data_filename,
        verbose=1):

    myVTK.myPrint(verbose, "*** readAbaqusStressFromDAT: " + data_filename + " ***")

    s_array = myVTK.createFloatArray("", 6)

    data_file = open(data_filename, 'r')
    context = ""
    k_cell = 0
    for line in data_file:
        if (context == "reading stresses"):
            #print line
            if ("MAXIMUM" in line):
                context = ""
                continue
            if ("OR" in line):
                splitted_line = line.split()
                assert (int(splitted_line[0]) == k_cell+1), "Wrong element number. Aborting."
                s_list = [float(splitted_line[k]) for k in xrange(3,9)]
                s_array.InsertNextTuple(s_list)
                k_cell += 1

        if (line == "    ELEMENT  PT FOOT-       S11         S22         S33         S12         S13         S23     \n"):
            context = "reading stresses"

    data_file.close()

    myVTK.myPrint(verbose, "n_tuples = " + str(s_array.GetNumberOfTuples()))

    return s_array





