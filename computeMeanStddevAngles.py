#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import math
import numpy
import vtk

import myVTKPythonLibrary as myVTK

########################################################################

def computeMeanStddevAngles(
        angles,
        angles_in_degrees=True,
        angles_in_pm_pi=True,
        verbose=1):

    myVTK.myPrint(verbose, "*** computeMeanStddevAngles ***")

    if (angles_in_degrees):
        if (angles_in_pm_pi):
            mean = math.atan2(numpy.mean([numpy.sin(2*numpy.array(angles)*numpy.pi/180)]),
                              numpy.mean([numpy.cos(2*numpy.array(angles)*numpy.pi/180)]))*180/math.pi/2
        else:
            mean = math.atan2(numpy.mean([numpy.sin(numpy.array(angles)*numpy.pi/180)]),
                              numpy.mean([numpy.cos(numpy.array(angles)*numpy.pi/180)]))*180/math.pi

        stddev = numpy.sqrt(numpy.mean(((((numpy.array(angles)-mean)+90)%180)-90)**2))
    else:
        if (angles_in_pm_pi):
            mean = math.atan2(numpy.mean([numpy.sin(2*numpy.array(angles))]),
                              numpy.mean([numpy.cos(2*numpy.array(angles))]))/2
        else:
            mean = math.atan2(numpy.mean([numpy.sin(numpy.array(angles))]),
                              numpy.mean([numpy.cos(numpy.array(angles))]))

        stddev = numpy.sqrt(numpy.mean(((((numpy.array(angles)-mean)+math.pi/2)%math.pi)-math.pi/2)**2))

    return (mean, stddev)

#def cleanAngles(
        #angles,
        #verbose=1):

    #myVTK.myPrint(verbose, "*** cleanAngles ***")

    #switch = True
    #while (switch):

        #avg = numpy.mean(angles)
        #if (avg < -90):
            #angles += 180
            #avg += 180
        #elif (avg > 90):
            #angles -= 180
            #avg -= 180

        #switch = False
        #for k in xrange(len(angles)):
            #if (angles[k] < avg-90.):
                #switch = True
                #angles[k] += 180.
            #if (angles[k] > avg+90.):
                #switch = True
                #angles[k] -= 180.
