#coding=utf8

########################################################################
###                                                                  ###
### Created by Martin Genet, 2012-2015                               ###
###                                                                  ###
### University of California at San Francisco (UCSF), USA            ###
### Swiss Federal Institute of Technology (ETH), Zurich, Switzerland ###
###                                                                  ###
########################################################################

import os

import myVTKPythonLibrary as myVTK

########################################################################

def exportDynaDeformationGradients(d3plot_file_basename):
    lspp_filename = d3plot_file_basename+'.lspp'
    lspp_file = open(lspp_filename, 'w')
    lspp_file.write('''\
open d3plot "''' + d3plot_file_basename + '''.d3plot"
selectpart shell off
fringe 91
output "''' + d3plot_file_basename + '''.history#11" 102 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 92
output "''' + d3plot_file_basename + '''.history#12" 102 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 93
output "''' + d3plot_file_basename + '''.history#13" 102 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 94
output "''' + d3plot_file_basename + '''.history#14" 102 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 95
output "''' + d3plot_file_basename + '''.history#15" 102 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 96
output "''' + d3plot_file_basename + '''.history#16" 102 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 97
output "''' + d3plot_file_basename + '''.history#17" 102 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 98
output "''' + d3plot_file_basename + '''.history#18" 102 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
fringe 99
output "''' + d3plot_file_basename + '''.history#19" 102 1 0 1 0 0 0 0 1 0 0 0 0 0 0 1.000000
exit
''')
    lspp_file.close()
    os.system('lspp -nographics c='+lspp_filename)
