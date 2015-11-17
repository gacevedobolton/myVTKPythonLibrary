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

import math
import random

import myVTKPythonLibrary as myVTK

########################################################################

def computeSyntheticHelixAngles2(
        farray_rr,
        farray_cc,
        farray_ll,
        angles_end=[[+60.], [+60.]],
        angles_epi=[[-60.], [-60.]],
        sigma=0.,
        farray_angle_helix=None,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeSyntheticHelixAngles2 ***")

    n_l = len(angles_end)
    assert (n_l > 1),\
        "n_l must be greater than 1. Aborting."
    assert (len(angles_epi) == n_l),\
        "angles_end and angle_epi must have same length (n_l). Aborting."
    d_l = 1./(n_l-1)

    n_c = len(angles_end[0])
    assert (n_c > 0),\
        "n_c must be greater than 0. Aborting."
    for angles in angles_end+angles_epi:
        assert (len(angles) == n_c),\
            "angles lists must have same length (n_c). Aborting."
    d_c = 1./n_c

    n_tuples = farray_rr.GetNumberOfTuples()
    assert (farray_cc.GetNumberOfTuples() == n_tuples)
    assert (farray_ll.GetNumberOfTuples() == n_tuples)

    if (farray_angle_helix == None):
        farray_angle_helix = myVTK.createFloatArray(
            "angle_helix",
            1,
            n_tuples)
    else:
        assert (farray_angle_helix.GetNumberOfTuples() == n_tuples)

    for k_tuple in xrange(n_tuples):
        #print "k_tuple = " + str(k_tuple)

        cc = farray_cc.GetTuple(k_tuple)[0]
        i_c = int(cc/d_c/1.000001)
        #print "i_c = " + str(i_c)

        zeta = (t - i_c*d_c) / d_c
        #print "zeta = " + str(zeta)

        ll = farray_ll.GetTuple(k_tuple)[0]
        i_l = int(ll/d_l/1.000001)
        #print "i_l = " + str(i_l)

        eta = (ll - i_l*d_l) / d_l
        #print "eta = " + str(eta)

        t_ii_end = angles_end[i_l  ][ i_c   %n_c]
        t_ji_end = angles_end[i_l  ][(i_c+1)%n_c]
        t_ij_end = angles_end[i_l+1][ i_c   %n_c]
        t_jj_end = angles_end[i_l+1][(i_c+1)%n_c]
        t_ii_epi = angles_epi[i_l  ][ i_c   %n_c]
        t_ji_epi = angles_epi[i_l  ][(i_c+1)%n_c]
        t_ij_epi = angles_epi[i_l+1][ i_c   %n_c]
        t_jj_epi = angles_epi[i_l+1][(i_c+1)%n_c]
        #print "t_ii_end = " + str(t_ii_end)
        #print "t_ji_end = " + str(t_ji_end)
        #print "t_ij_end = " + str(t_ij_end)
        #print "t_jj_end = " + str(t_jj_end)
        #print "t_ii_epi = " + str(t_ii_epi)
        #print "t_ji_epi = " + str(t_ji_epi)
        #print "t_ij_epi = " + str(t_ij_epi)
        #print "t_jj_epi = " + str(t_jj_epi)

        helix_angle_end = t_ii_end * (1 - zeta - eta + zeta*eta) \
                        + t_ji_end * (zeta - zeta*eta) \
                        + t_ij_end * (eta - zeta*eta) \
                        + t_jj_end * (zeta*eta)
        helix_angle_epi = t_ii_epi * (1 - zeta - eta + zeta*eta) \
                        + t_ji_epi * (zeta - zeta*eta) \
                        + t_ij_epi * (eta - zeta*eta) \
                        + t_jj_epi * (zeta*eta)

        rr = farray_rr.GetTuple(k_tuple)[0]
        helix_angle_in_degrees = (1.-rr) * helix_angle_end \
                               +     rr  * helix_angle_epi

        if (sigma > 0.):
            helix_angle_in_degrees += random.normalvariate(0., sigma)
            helix_angle_in_degrees  = (helix_angle_in_degrees+90.)%180.-90.

        farray_angle_helix.InsertTuple(
            k_tuple,
            [helix_angle_in_degrees])

    return farray_angle_helix

########################################################################

def addSyntheticHelixAngles2(
        ugrid,
        angles_end,
        angles_epi,
        type_of_support="cell",
        sigma=0,
        verbose=1):

    myVTK.myPrint(verbose, "*** addSyntheticHelixAngles2 ***")

    if (type_of_support == "cell"):
        ugrid_data = ugrid.GetCellData()
    elif (type_of_support == "point"):
        ugrid_data = ugrid.GetPointData()

    farray_rr = ugrid_data.GetArray("rr")
    farray_cc = ugrid_data.GetArray("cc")
    farray_ll = ugrid_data.GetArray("ll")

    farray_angle_helix = ugrid_data.GetArray("angle_helix")

    farray_angle_helix = computeSyntheticHelixAngles2(
        farray_rr=farray_rr,
        farray_cc=farray_cc,
        farray_ll=farray_ll,
        angles_end=angles_end,
        angles_epi=angles_epi,
        sigma=sigma,
        farray_angle_helix=farray_angle_helix,
        verbose=verbose-1)

    ugrid_data.AddArray(farray_angle_helix)

    return farray_angle_helix
