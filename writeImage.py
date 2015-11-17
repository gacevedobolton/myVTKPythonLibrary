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

import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def writeImage(
        image,
        filename,
        verbose=1):

    myVTK.myPrint(verbose, "*** writeImage: " + filename + " ***")

    if ('vtk' in filename):
        image_writer = vtk.vtkImageWriter()
    elif ('vti' in filename):
        image_writer = vtk.vtkXMLImageWriter()
    else:
        assert 0, "File must be .vtk or .vti. Aborting."

    image_writer.SetFileName(filename)
    image_writer.SetInputData(image)
    image_writer.Update()
    image_writer.Write()
