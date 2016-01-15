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

import myVTKPythonLibrary as myVTK

########################################################################

def getPDataSurfaceArea(
        pdata,
        verbose=1):

    myVTK.myPrint(verbose, "*** getPDataSurfaceArea ***")

    mass_properties = myVTK.createMassProperties(
        pdata=pdata,
        verbose=verbose-1)

    return mass_properties.GetSurfaceArea()
