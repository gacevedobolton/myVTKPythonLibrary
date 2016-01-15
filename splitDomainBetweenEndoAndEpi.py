#!/usr/bin/python
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

import argparse
import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def splitDomainBetweenEndoAndEpi(
        pdata_domain,
        ratio=0.90,
        verbose=1):

    myVTK.myPrint(verbose, "*** splitDomainBetweenEndoAndEpi ***")

    bounds = pdata_domain.GetBounds()
    C = [(bounds[0]+bounds[1])/2, (bounds[2]+bounds[3])/2, (1.-ratio)*bounds[4]+ratio*bounds[5]]
    N = [0,0,1]

    (pdata_domain,
     cap) = myVTK.clipPDataUsingPlane(
         pdata_mesh=pdata_domain,
         plane_C=C,
         plane_N=N,
         verbose=verbose-1)

    connectivity0 = vtk.vtkPolyDataConnectivityFilter()
    connectivity0.SetExtractionModeToSpecifiedRegions()
    connectivity0.AddSpecifiedRegion(0)
    connectivity0.SetInputData(pdata_domain)
    connectivity0.Update()
    pdata0 = connectivity0.GetOutput()

    connectivity1 = vtk.vtkPolyDataConnectivityFilter()
    connectivity1.SetExtractionModeToSpecifiedRegions()
    connectivity1.AddSpecifiedRegion(1)
    connectivity1.SetInputData(pdata_domain)
    connectivity1.Update()
    pdata1 = connectivity1.GetOutput()

    if (myVTK.getPDataSurfaceArea(pdata0,0) < myVTK.getPDataSurfaceArea(pdata1,0)):
        return pdata0, pdata1
    else:
        return pdata1, pdata0

########################################################################

if (__name__ == "__main__"):
    parser = argparse.ArgumentParser()
    parser.add_argument('domain_filename' , type=str              )
    parser.add_argument('--endLV_filename', type=str, default=None)
    parser.add_argument('--epiLV_filename', type=str, default=None)
    parser.add_argument('--verbose', '-v' , type=int, default=1   )
    args = parser.parse_args()

    if (args.endLV_filename == None):
        args.endLV_filename = args.domain_filename.replace("LV", "EndLV")
    if (args.epiLV_filename == None):
        args.epiLV_filename = args.domain_filename.replace("LV", "EpiLV")

    pdata_domain = myVTK.readSTL(
        filename=args.domain_filename,
        verbose=args.verbose)

    (pdata_end,
     pdata_epi) = myVTK.splitDomainBetweenEndoAndEpi(
         pdata_domain=pdata_domain,
         verbose=args.verbose)

    myVTK.writeSTL(
        pdata=pdata_end,
        filename=args.endLV_filename,
        verbose=args.verbose)

    myVTK.writeSTL(
        pdata=pdata_epi,
        filename=args.epiLV_filename,
        verbose=args.verbose)
