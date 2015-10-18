#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import sys
import math
import random
import numpy

import myVTKPythonLibrary as myVTK

########################################################################

def writeFiberOrientationFileForGNUPlot(
        angles_end,
        angles_epi,
        filename,
        verbose=1):

    myVTK.myPrint(verbose, "*** writeFiberOrientationFileForGNUPlot ***")

    assert (len(angles_end) == len(angles_epi)), "angles_end and angle_epi must have same length (n_l). Aborting."
    n_l = len(angles_end)
    d_l = 1./(n_l-1)
    n_c = len(angles_end[0])
    for angles in angles_end+angles_epi:
        assert (len(angles) == n_c), "angles lists must have same length (n_c). Aborting."
    d_c = 1./n_c

    fiber_orientation_file = open(filename, 'w')
    fiber_orientation_file.write('# c ang_end ang_epi z\n')

    n_c = 12
    for k_c in xrange(n_c+1):
        c    = float(k_c) / n_c
        i_c  = int(c/d_c/1.000001)
        zeta = (c - i_c*d_c) / d_c

        n_l = 10
        for k_l in xrange(n_l+1):
            l   = float(k_l) / n_l
            i_l = int(l/d_l/1.000001)
            eta = (l - i_l*d_l) / d_l

            t_ii_end = angles_end[i_l  ][ i_c   %n_c]
            t_ji_end = angles_end[i_l  ][(i_c+1)%n_c]
            t_ij_end = angles_end[i_l+1][ i_c   %n_c]
            t_jj_end = angles_end[i_l+1][(i_c+1)%n_c]
            t_ii_epi = angles_epi[i_l  ][ i_c   %n_c]
            t_ji_epi = angles_epi[i_l  ][(i_c+1)%n_c]
            t_ij_epi = angles_epi[i_l+1][ i_c   %n_c]
            t_jj_epi = angles_epi[i_l+1][(i_c+1)%n_c]

            helix_angle_end = t_ii_end * (1 - zeta - eta + zeta*eta) \
                            + t_ji_end * (    zeta       - zeta*eta) \
                            + t_ij_end * (           eta - zeta*eta) \
                            + t_jj_end * (                 zeta*eta)
            helix_angle_epi = t_ii_epi * (1 - zeta - eta + zeta*eta) \
                            + t_ji_epi * (    zeta       - zeta*eta) \
                            + t_ij_epi * (           eta - zeta*eta) \
                            + t_jj_epi * (                 zeta*eta)

            fiber_orientation_file.write(" ".join([str(x) for x in [t, helix_angle_end, helix_angle_epi, z]]) + "\n")

        fiber_orientation_file.write("\n")

    fiber_orientation_file.close()
