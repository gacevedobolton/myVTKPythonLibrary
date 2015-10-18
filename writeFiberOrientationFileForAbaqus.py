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

def writeFiberOrientationFileForAbaqus(
        mesh,
        filename,
        eF_field_name="eF",
        eS_field_name="eS",
        sep=", ",
        verbose=1):

    myVTK.myPrint(verbose, "*** writeFiberOrientationFileForAbaqus ***")

    orientation_file = open(filename, "w")
    orientation_file.write(", 1., 0., 0., 0., 1., 0." + "\n")

    n_cells = mesh.GetNumberOfCells()

    eF_array = mesh.GetCellData().GetArray(eF_field_name)
    eS_array = mesh.GetCellData().GetArray(eS_field_name)

    for k_cell in xrange(n_cells):
        eF = eF_array.GetTuple(k_cell)
        eS = eS_array.GetTuple(k_cell)

        line = str(k_cell+1)
        for k in xrange(3): line += sep + str(eF[k])
        for k in xrange(3): line += sep + str(eS[k])
        line += "\n"
        orientation_file.write(line)

    orientation_file.close()
