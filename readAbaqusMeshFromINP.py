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

import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def readAbaqusMeshFromINP(
        mesh_filename,
        elem_types="all",
        verbose=1):

    myVTK.myPrint(verbose, "*** readAbaqusMeshFromINP: " + mesh_filename + " ***")

    points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()

    mesh_file = open(mesh_filename, "r")

    context = ""
    for line in mesh_file:
        if (line[-1:] == "\n"): line = line[:-1]
        #myVTK.myPrint(verbose, "line =", line)

        if line.startswith("**"): continue

        if (context == "reading nodes"):
            if line.startswith("*"):
                context = ""
            else:
                splitted_line = line.split(",")
                points.InsertNextPoint([float(coord) for coord in splitted_line[1:4]])

        if (context == "reading elems"):
            if line.startswith("*"):
                context = ""
            else:
                splitted_line = line.split(",")
                assert (len(splitted_line) == 1+cell_n_points), "Wrong number of elements in line. Aborting."
                for k_point in xrange(cell_n_points): cell.GetPointIds().SetId(k_point, int(splitted_line[1+k_point])-1)
                cell_array.InsertNextCell(cell)

        if line.startswith("*NODE"):
            context = "reading nodes"
        if line.startswith("*ELEMENT"):
            if ("TYPE=F3D4" in line) and (("quad" in elem_types) or ("all" in elem_types)):
                context = "reading elems"
                cell_vtk_type = vtk.VTK_QUAD
                cell_n_points = 4
                cell = vtk.vtkQuad()
            elif ("TYPE=C3D4" in line) and (("tet" in elem_types) or ("all" in elem_types)):
                context = "reading elems"
                cell_vtk_type = vtk.VTK_TETRA
                cell_n_points = 4
                cell = vtk.vtkTetra()
            elif ("TYPE=C3D8" in line) and (("hex" in elem_types) or ("all" in elem_types)):
                context = "reading elems"
                cell_vtk_type = vtk.VTK_HEXAHEDRON
                cell_n_points = 8
                cell = vtk.vtkHexahedron()
            else:
                print "Warning: elements not read: " + line + "."

    mesh_file.close()

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.SetCells(cell_vtk_type, cell_array)

    myVTK.myPrint(verbose, "n_cells = " + str(ugrid.GetNumberOfCells()))

    return ugrid
