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

def readDynaMesh(lsdyna_mesh_filename,
                 cell_type='hexahedron',
                 verbose=1):

    myVTK.myPrint(verbose, "*** readDynaMesh: " + lsdyna_mesh_filename + " ***")

    points = vtk.vtkPoints()

    if (cell_type == 'vertex'):
        cell_vtk_type = vtk.VTK_VERTEX
        cell = vtk.vtkVertex()
        cell_array = vtk.vtkCellArray()
    elif (cell_type == 'hexahedron'):
        cell_vtk_type = vtk.VTK_HEXAHEDRON
        cell = vtk.vtkHexahedron()
        cell_array = vtk.vtkCellArray()
    else:
        print 'Wrong cell type. Aborting.'
        exit()

    myVTK.myPrint(verbose, "Reading Dyna mesh file...")

    mesh_file = open(lsdyna_mesh_filename, 'r')

    context = ''
    for line in mesh_file:
        if (line[-1:] == '\n'): line = line[:-1]
        #myVTK.myPrint(verbose, "line ="+line)

        if line.startswith('$'): continue

        if (context == 'reading nodes'):
            if line.startswith('*'):
                context = ''
            else:
                splitted_line = line.split(',')
                points.InsertNextPoint([float(coord) for coord in splitted_line[1:4]])
                if (cell_type == 'vertex'):
                    cell.GetPointIds().SetId(0, points.GetNumberOfPoints()-1)
                    cell_array.InsertNextCell(cell)

        if (context == 'reading elems'):
            if line.startswith('*'):
                context = ''
            else:
                splitted_line = line.split(',')
                if (len(splitted_line) == 3): continue
                if (cell_type == 'hexahedron'):
                    for k_node in xrange(8):
                        cell.GetPointIds().SetId(k_node, int(splitted_line[2+k_node])-1)
                cell_array.InsertNextCell(cell)

        if line.startswith('*NODE'): context = 'reading nodes'
        if line.startswith('*ELEMENT_SOLID'): context = 'reading elems'

    mesh_file.close()

    myVTK.myPrint(verbose, "Creating VTK mesh...")

    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(points)
    ugrid.SetCells(cell_vtk_type, cell_array)

    return ugrid
