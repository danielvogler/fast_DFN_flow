#! /usr/bin/python

# Copyright:
# 	Alex Hobe
# 	Daniel Vogler
# 	Martin P. Seybold




#########
### This script runs the entire developed framework.
#########


import time

saving 				= 0		# save results to file

segment 			= 0 	# use HSPM
simpleWidth 		= 1     # use ISPM

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
