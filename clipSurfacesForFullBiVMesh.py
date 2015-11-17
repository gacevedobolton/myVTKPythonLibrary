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

import argparse
import sys
import vtk

import myVTKPythonLibrary as myVTK
from mat_vec_tools import *

########################################################################

def clipSurfacesForFullBiVMesh(
        pdata_endLV,
        pdata_endRV,
        pdata_epi,
        verbose=1):

    myVTK.myPrint(verbose, "*** clipSurfacesForFullBiVMesh ***")

    pdata_endLV_implicit_distance = vtk.vtkImplicitPolyDataDistance()
    pdata_endLV_implicit_distance.SetInput(pdata_endLV)

    pdata_endRV_implicit_distance = vtk.vtkImplicitPolyDataDistance()
    pdata_endRV_implicit_distance.SetInput(pdata_endRV)

    pdata_epi_implicit_distance = vtk.vtkImplicitPolyDataDistance()
    pdata_epi_implicit_distance.SetInput(pdata_epi)

    pdata_endLV_clip = vtk.vtkClipPolyData()
    pdata_endLV_clip.SetInputData(pdata_endLV)
    pdata_endLV_clip.SetClipFunction(pdata_epi_implicit_distance)
    pdata_endLV_clip.InsideOutOn()
    pdata_endLV_clip.Update()
    clipped_pdata_endLV = pdata_endLV_clip.GetOutput(0)

    pdata_endRV_clip = vtk.vtkClipPolyData()
    pdata_endRV_clip.SetInputData(pdata_endRV)
    pdata_endRV_clip.SetClipFunction(pdata_epi_implicit_distance)
    pdata_endRV_clip.InsideOutOn()
    pdata_endRV_clip.Update()
    clipped_pdata_endRV = pdata_endRV_clip.GetOutput(0)

    pdata_epi_clip = vtk.vtkClipPolyData()
    pdata_epi_clip.SetInputData(pdata_epi)
    pdata_epi_clip.SetClipFunction(pdata_endLV_implicit_distance)
    pdata_epi_clip.GenerateClippedOutputOn()
    pdata_epi_clip.Update()
    clipped_pdata_epi = pdata_epi_clip.GetOutput(0)
    clipped_valM = pdata_epi_clip.GetOutput(1)

    pdata_epi_clip = vtk.vtkClipPolyData()
    pdata_epi_clip.SetInputData(clipped_pdata_epi)
    pdata_epi_clip.SetClipFunction(pdata_endRV_implicit_distance)
    pdata_epi_clip.GenerateClippedOutputOn()
    pdata_epi_clip.Update()
    clipped_pdata_epi = pdata_epi_clip.GetOutput(0)
    clipped_valP = pdata_epi_clip.GetOutput(1)

    return (clipped_pdata_endLV,
            clipped_pdata_endRV,
            clipped_pdata_epi,
            clipped_valM,
            clipped_valP)

########################################################################

if (__name__ == "__main__"):

    parser = argparse.ArgumentParser()
    parser.add_argument('endLV_filename', type=str)
    parser.add_argument('endRV_filename', type=str)
    parser.add_argument('epi_filename', type=str)
    parser.add_argument('-v', '--verbose', type=int, default=1)
    args = parser.parse_args()

    pdata_endLV = myVTK.readSTL(
        filename=args.endLV_filename,
        verbose=args.verbose)
    pdata_endRV = myVTK.readSTL(
        filename=endRV_filename,
        verbose=args.verbose)
    pdata_epi = myVTK.readSTL(
        filename=args.epi_filename,
        verbose=args.verbose)

    (clipped_pdata_endLV,
     clipped_pdata_endRV,
     clipped_pdata_epi,
     clipped_valM,
     clipped_valP) = myVTK.clipSurfacesForFullBiVMesh(
        pdata_endLV=pdata_endLV,
        pdata_endRV=pdata_endRV,
        pdata_epi=pdata_epi,
        verbose=args.verbose)

    myVTK.writeSTL(
        pdata=clipped_pdata_endLV,
        filename="endLV.stl",
        verbose=args.verbose)
    myVTK.writeSTL(
        pdata=clipped_pdata_endRV,
        filename="endRV.stl",
        verbose=args.verbose)
    myVTK.writeSTL(
        pdata=clipped_pdata_epi,
        filename="epi.stl",
        verbose=args.verbose)
    myVTK.writeSTL(
        pdata=clipped_pdata_valM,
        filename="valM.stl",
        verbose=args.verbose)
    myVTK.writeSTL(
        pdata=clipped_pdata_valP,
        filename="valP.stl",
        verbose=args.verbose)

