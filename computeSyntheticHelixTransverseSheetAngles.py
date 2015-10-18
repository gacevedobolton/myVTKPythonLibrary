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

def computeSyntheticHelixTransverseSheetAngles(
        farray_rr,
        farray_cc,
        farray_ll,
        angles="+/-60",
        sigma=0.,
        farray_angle_helix=None,
        farray_angle_trans=None,
        farray_angle_sheet=None,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeSyntheticHelixTransverseSheetAngles ***")

    if (angles == "+/-60"):
        angles = [[[[+60], [+60]], [[-60], [-60]]],
                  [[[  0], [  0]], [[  0], [  0]]],
                  [[[  0], [  0]], [[  0], [  0]]]]

    n_l = [[0. for k_r in xrange(2)] for k_angle in xrange(3)]
    d_l = [[0. for k_r in xrange(2)] for k_angle in xrange(3)]
    n_c = [[0. for k_r in xrange(2)] for k_angle in xrange(3)]
    d_c = [[0. for k_r in xrange(2)] for k_angle in xrange(3)]
    for k_angle in xrange(3):
        #print "k_angle = " + str(k_angle)
        for k_r in xrange(2):
            #print "k_r = " + str(k_r)

            n_l[k_angle][k_r] = len(angles[k_angle][k_r])
            assert (n_l[k_angle][k_r] > 1), "Must have more than 1 longitudinal dof. Aborting."
            #print "n_l = " + str(n_l)
            d_l[k_angle][k_r] = 1./(n_l[k_angle][k_r]-1)
            #print "d_l = " + str(d_l)

            n_c[k_angle][k_r] = len(angles[k_angle][k_r][0])
            assert (n_c[k_angle][k_r] > 0), "Must have more than 0 circumferential dof. Aborting."
            #print "n_c = " + str(n_c)
            d_c[k_angle][k_r] = 1./n_c[k_angle][k_r]
            #print "d_c = " + str(d_c)

    n_tuples = farray_rr.GetNumberOfTuples()

    if (farray_angle_helix == None):
        farray_angle_helix = myVTK.createFloatArray(
            "angle_helix",
            1,
            n_tuples)
    if (farray_angle_trans == None):
        farray_angle_trans = myVTK.createFloatArray(
            "angle_trans",
            1,
            n_tuples)
    if (farray_angle_sheet == None):
        farray_angle_sheet = myVTK.createFloatArray(
            "angle_sheet",
            1,
            n_tuples)
    farray_angles = [farray_angle_helix,\
                     farray_angle_trans,\
                     farray_angle_sheet]

    for k_tuple in xrange(n_tuples):
        #print "k_tuple = " + str(k_tuple)

        angles_in_degrees = numpy.array([0.]*3)
        for k_angle in xrange(3):
            #print "k_angle = " + str(k_angle)

                rr = farray_rr.GetTuple(k_tuple)[0]

                cc = farray_cc.GetTuple(k_tuple)[0]
                i_c_end = int(cc/d_c[k_angle][0]/1.000001)
                i_c_epi = int(cc/d_c[k_angle][1]/1.000001)
                #print "i_c_end = " + str(i_c_end)
                #print "i_c_epi = " + str(i_c_epi)

                zeta_end = (cc - i_c_end*d_c[k_angle][0])/d_c[k_angle][0]
                zeta_epi = (cc - i_c_epi*d_c[k_angle][1])/d_c[k_angle][1]
                #print "zeta_end = " + str(zeta_end)
                #print "zeta_epi = " + str(zeta_epi)

                ll = farray_ll.GetTuple(k_tuple)[0]
                i_l_end = int(ll/d_l[k_angle][0]/1.000001)
                i_l_epi = int(ll/d_l[k_angle][1]/1.000001)
                #print "i_l_end = " + str(i_l_end)
                #print "i_l_epi = " + str(i_l_epi)

                eta_end = (ll - i_l_end*d_l[k_angle][0])/d_l[k_angle][0]
                eta_epi = (ll - i_l_epi*d_l[k_angle][1])/d_l[k_angle][1]
                #print "eta_end = " + str(eta_end)
                #print "eta_epi = " + str(eta_epi)

                t_ii_end = angles[k_angle][0][i_l_end  ][ i_c_end   %n_c[k_angle][0]]
                t_ji_end = angles[k_angle][0][i_l_end  ][(i_c_end+1)%n_c[k_angle][0]]
                t_ij_end = angles[k_angle][0][i_l_end+1][ i_c_end   %n_c[k_angle][0]]
                t_jj_end = angles[k_angle][0][i_l_end+1][(i_c_end+1)%n_c[k_angle][0]]
                t_ii_epi = angles[k_angle][1][i_l_epi  ][ i_c_epi   %n_c[k_angle][1]]
                t_ji_epi = angles[k_angle][1][i_l_epi  ][(i_c_epi+1)%n_c[k_angle][1]]
                t_ij_epi = angles[k_angle][1][i_l_epi+1][ i_c_epi   %n_c[k_angle][1]]
                t_jj_epi = angles[k_angle][1][i_l_epi+1][(i_c_epi+1)%n_c[k_angle][1]]
                #print "t_ii_end = " + str(t_ii_end)
                #print "t_ji_end = " + str(t_ji_end)
                #print "t_ij_end = " + str(t_ij_end)
                #print "t_jj_end = " + str(t_jj_end)
                #print "t_ii_epi = " + str(t_ii_epi)
                #print "t_ji_epi = " + str(t_ji_epi)
                #print "t_ij_epi = " + str(t_ij_epi)
                #print "t_jj_epi = " + str(t_jj_epi)

                angle_end = t_ii_end * (1 - zeta_end - eta_end + zeta_end*eta_end) \
                          + t_ji_end * (    zeta_end           - zeta_end*eta_end) \
                          + t_ij_end * (               eta_end - zeta_end*eta_end) \
                          + t_jj_end * (                         zeta_end*eta_end)
                angle_epi = t_ii_epi * (1 - zeta_epi - eta_epi + zeta_epi*eta_epi) \
                          + t_ji_epi * (    zeta_epi           - zeta_epi*eta_epi) \
                          + t_ij_epi * (               eta_epi - zeta_epi*eta_epi) \
                          + t_jj_epi * (                         zeta_epi*eta_epi)

                angles_in_degrees[k_angle] = (1.-rr) * angle_end \
                                           +     rr  * angle_epi
                if (sigma > 0.):
                    angles_in_degrees[k_angle] += random.normalvariate(0., sigma)
                    angles_in_degrees[k_angle]  = (angles_in_degrees[k_angle]+90)%180-90
                farray_angles[k_angle].InsertTuple(k_tuple, [angles_in_degrees[k_angle]])

    return (farray_angle_helix,
            farray_angle_trans,
            farray_angle_sheet)

########################################################################

def addSyntheticHelixTransverseSheetAngles(
        ugrid,
        angles="+/-60",
        type_of_support="cell",
        sigma=0.,
        verbose=1):

    myVTK.myPrint(verbose, "*** addSyntheticHelixTransverseSheetAngles ***")

    if (type_of_support == "cell"):
        ugrid_data = ugrid.GetCellData()
    elif (type_of_support == "point"):
        ugrid_data = ugrid.GetPointData()

    farray_rr = ugrid_data.GetArray("rr")
    farray_cc = ugrid_data.GetArray("cc")
    farray_ll = ugrid_data.GetArray("ll")

    farray_angle_helix = ugrid_data.GetArray("angle_helix")
    farray_angle_trans = ugrid_data.GetArray("angle_trans")
    farray_angle_sheet = ugrid_data.GetArray("angle_sheet")

    (farray_angle_helix,
     farray_angle_trans,
     farray_angle_sheet) = computeSyntheticHelixTransverseSheetAngles(
        farray_rr=farray_rr,
        farray_cc=farray_cc,
        farray_ll=farray_ll,
        angles=angles,
        sigma=sigma,
        farray_angle_helix=farray_angle_helix,
        farray_angle_trans=farray_angle_trans,
        farray_angle_sheet=farray_angle_sheet,
        verbose=verbose-1)

    ugrid_data.AddArray(farray_angle_helix)
    ugrid_data.AddArray(farray_angle_trans)
    ugrid_data.AddArray(farray_angle_sheet)

    return (farray_angle_helix,
            farray_angle_trans,
            farray_angle_sheet)
