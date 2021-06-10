#! /usr/bin/python

# Copyright:
# 	Alex Hobe
# 	Daniel Vogler
# 	Martin P. Seybold


# TODO
# - imports
# - check out inheritence

def visualizeGraphPipeline( graphFileName ):
	# plotting of 2D and 3D examples is equivalent
	g = loadGraphWithValues( graphFileName )
	flow = calculateFlow( g )
	plotGraph( g )
	plotMinCutFlow( g, flow )		
	plotPathsOnGraph( g, flow )

def plotPathsOnGraph( g, flow ):
	# This one overlaps greatly with the flow-rate estimation

	# while shortestPath exists
		nextPath = getPath( g, flow )
		plotPath( g, nextPath )
		removePath( g, nextpath )
