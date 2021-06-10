############################
#
# Copyright (c) 2015, 2016, 2017
# ETH Swiss Federal Institute of Technology Zurich, Zurich, Switzerland and
# University of Stuttgart, Stuttgart, Germany
#
# All rights reserved.
#
# Fast Methods in Discrete Fracture Networks, Version 1.0.0
#
# Written by:
#     Martin P. Seybold - seybold(at)fmi.uni-stuttgart.de
#     Daniel Vogler - davogler(at)ethz.ch
#     Alex Hobe - ahobe(at)student.ethz.ch
#
# The copyright to the software herein is the property of
# Daniel Vogler and Alex Hobe of ETH Swiss Federal Institute
# of Technology Zurich, Zurich, Switzerland and Martin P. Seybold
# of University of Stuttgart, Stuttgart, Germany.
# The software may be used only with the written permission
# of Martin P. Seybold or Daniel Vogler or in accordance with
# the terms and conditions stipulated in the agreement/contract
# under which the software has been supplied.
#
# Software is provided as is, without warranties of any kind,
# express or implied, including warranties of merchantability or
# fitness for a particular purpose. Software is provided without
# representations of authorship or originality.
#
# This copyright notice must not be removed.
#
############################

#! /usr/bin/python

__author__ = "martin"
__date__ = "$Jun 21, 2015 11:00:39 PM$"


#########
### This script runs the entire developed framework.
#########


import time

saving 				= 0		# save results to file

segment 			= 1 	# use HSPM
simpleWidth 		= 0     # use ISPM

if saving == 1:
	filename = 'ISPM_100_case_results_Quartz.txt'
	#filename = 'HSPM_100_case_results_Quartz.txt'
	fileID = open(filename,'w'); # open file for writing; overwrite existing contents
	fileID.write("simulation, Qgrav, numpaths, runtime [s], nVertices, nEdges  \n")

for i in range(66,69):
	##########################
	### Choosing cases
	##########################
	if i > 65  and i < 68:
		workingDir = "../benchmark/bm54_f17/plot/"
		inputFile  =              "bm54_f17_3.xml"
	if i > 67 and i < 69:
		workingDir = "../benchmark/bm95_f250/plot/"
		inputFile  =              "bm95_f250_0_2.xml"

	############################
	### Running the algorithm
	############################

	t = time.time()



	if segment == 1:
		from HSPM_opt import *
		Q, numPaths, nVertices, nEdges, tRead, tG1, tG2, tG3, tG4, tGraph, tMaxFlow, tBugFix, tPath  = segmentation(workingDir, inputFile)

	if simpleWidth == 1:
		from ISPM_opt import *
		Q, numPaths, nVertices, nEdges, tRead, tG1, tG2, tG3, tG4, tGraph, tMaxFlow, tBugFix, tPath = simpleMethod(workingDir, inputFile)

	elapsed = time.time() - t

	print("Input file: {}".format( inputFile ) )


	if saving == 1 and i>66:
			fileID.write(inputFile +  ","  +"{0},".format(Q) + "{0},".format(numPaths) + "{0:.2f},".format(elapsed) + "{0},".format(nVertices) + "{0} \n".format(nEdges) )


	print("Flow rate: {}".format( Q ) )
	print("Number of paths: {}".format( numPaths ) )
	print("{0:.2f}".format(elapsed) )

	if saving == 1:
		fileID.close()
