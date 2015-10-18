#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def readImage(
        filename,
        verbose=1):

    myVTK.myPrint(verbose, "*** readImage: " + filename + " ***")

    if ('vtk' in filename):
        image_reader = vtk.vtkImageDataReader()
    elif ('vti' in filename):
        image_reader = vtk.vtkXMLImageDataReader()
    else:
        assert 0, "File must be .vtk or .vti. Aborting."

    image_reader.SetFileName(filename)
    image_reader.Update()
    image = image_reader.GetOutput()

    if (verbose):
        print "n_points = " + str(image.GetNumberOfPoints())

    return image
