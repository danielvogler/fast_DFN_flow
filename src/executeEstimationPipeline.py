#! /usr/bin/python

# Copyright:
# 	Alex Hobe
# 	Daniel Vogler
# 	Martin P. Seybold


# TODO
# - imports
# - check out inheritence

def executeEstimationPipeline(fileName,graphType):
	fileType  = getFileType(fileName)
	fracBoxes = getFractureBoxes(fileName,fileType)
	
	# Point of abstraction
	# Everything after this point assumes the same structure of fracBoxes
	
	fracGraph = makeGraph(fracBoxes, graphType)

	Q = estimateFlowRate(fracGraph)
	
	return Q
	

# This would allow an abstract pipeline class using the initial stages.
# filetype, fracBoxes, and fracGraph.
# The following classes can inherit from them
# estimationPipeline
# plottingPipeline
# saveGraphPipeline
# other?

# Saving a graph would allow the plotting pipeline to ignore the initial stage 
# and just use readGraph() and plotting routines.


# network creator
# network investigator

# test driven development?
