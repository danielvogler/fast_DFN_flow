#! /usr/bin/python

# Copyright:
# 	Alex Hobe
# 	Daniel Vogler
# 	Martin P. Seybold


# TODO
# - imports
# - refactor further

import sys
import numpy as np
from graphToolwrapper import *

gt = graphToolwrapper.all()
tp = graphToolwrapper.topology()


def abstractQestimator(fracGraph, solver):
	""" Set of methods, that estimate flowrates of a given graph.

	"""

def estimateFlowRate( g, solver ):
	# requires cap and eps
	cap = g.edge_properties["cap"]

	if solver == "HSPM": # shortest path maxflow		
		g, max_flow = solveMaxFlow(cap)

		# terminate of no percolation
		if abs(max_flow - 0.0)<eps:
			print( "" )
			print( "Graph not connected. No percolation possible." )
			print( "" )
			return 0.0, 0, nVertices, nEdges # Qvalue, number of paths.

		removeNoFlowEdges(g, res, eps)

		g, res = singleDirectionFlowCleanup( g, res )

		removeNoFlowEdges(g, res, eps)

		QgravTot, numPaths   = pathFinder( g, src, tgt, res, path_crit, aCmC, eps )
		return QgravTot, numPaths, nVertices, nEdges

	else:
		print( " flowrate solver: ' ' not yet implemented" )
		# TODO: error/exception instead of print
		sys.exit()

	return Q


def solveMaxFlow(g):
	cap           = g.edge_properties["cap"]  
	full_source = 0
	full_target = 1
	src, tgt = g.vertex(full_source), g.vertex(full_target)               # set the source and target nodes
	res = gt.edmonds_karp_max_flow(g, src, tgt, cap)
	res.a = cap.a - res.a  # the actual flow
	max_flow = sum(res[e] for e in tgt.in_edges())    # the maximum flow is the sum of the individual paths

	return g, max_flow, res, src, tgt 


def removeNoFlowEdges(g, res, eps):
	noFlow = gt.find_edge_range(g, res, [0.0, eps])
	for s in noFlow:
	    g.remove_edge(s)


def singleDirectionFlowCleanup( g, res ):
	# Clean up, so flow only goes in one direction
	for e in g.edges():
		eSrc = e.source()
		eTgt = e.target()
		if not g.edge(eTgt,eSrc):
			continue
		else:
			eReverse = g.edge(eTgt,eSrc)

		there = res[e]
		back  = res[eReverse]

		if there > 0 and back > 0:
			if there >= back:
				res[e] -= back
				res[eReverse] = 0
			else:
				res[eReverse] -= there
				res[e] = 0
	return g, res





def pathFinder( gCopy, src, tgt, resCopy, path_crit, aCmC, eps ):
	# TODO
	# - fix indentation
	# - test resCopy.copy()

	cent = gCopy.vertex_properties["cent"]
	e_length      = gCopy.edge_properties["e_length"]
	e_width       = gCopy.edge_properties["e_width"]
	# cap           = gCopy.edge_properties["cap"]      # capacity for flow in each edge
	# path_crit     = gCopy.edge_properties["path_crit"]

	numPaths, resFull, QgravTot, gravPart = initializePathFinder(resCopy)

	########
	# change this when using a different DFN size
	dxPart   = (10.6)**2
	########

	for i in range(0,1000): # perhaps do a while loop instead. This should be safer...

		vShort, eShort = tp.shortest_path(gCopy,src,tgt, weights=path_crit) # find shortest path

		if len(eShort) == 0: # no more path found
		    return QgravTot, numPaths

		else:
			pLength, pWidth, pFactor, pRes, minFlow, wSuperWeighted, Ltot = initializeIndividualPath()
			
			eStart, eEnd, eShort = getAndRemoveStartAndEndFromPath(eShort)

			minFlow, Ltot, wSuperWeighted = getMinFlowAndWeightsForHananPath( eShort, e_length, Ltot, wSuperWeighted, resFull, e_width, resCopy, minFlow )

			pGrav, Ldiag = getGravitationalAndDiagonalPathParts(vShort, cent, gravPart)

			resCopy = removeFlowFromPath( eStart, eEnd, eShort, resCopy, minFlow )

			QgravTot = addQforPath( aCmC, wSuperWeighted, Ldiag, pGrav, QgravTot )			

			removeNoFlowEdges(gCopy, resCopy, eps)

			numPaths += 1


def initializePathFinder(resCopy):
	numPaths = 0
	resFull  = resCopy.copy()
	QgravTot = 0.
	gravPart = 9.81*1000

	return numPaths, resFull, QgravTot, gravPart


def initializeIndividualPath():
	pLength = []  # store lengths of this path
	pWidth  = []  # store widths of this path
	pFactor = []  # store flow factor of this path (percentage of the edge that this path uses).
	pRes    = []  # store the flow results for this path
	minFlow = 10000
	wSuperWeighted = 0.0
	Ltot = 0.

	return pLength, pWidth, pFactor, pRes, minFlow, wSuperWeighted, Ltot


def getAndRemoveStartAndEndFromPath(eShort):
	eStart = eShort[0]
	eEnd   = eShort[-1]
	del eShort[0]
	del eShort[-1]

	return eStart, eEnd, eShort

def getMinFlowAndWeightsForHananPath( eShort, e_length, Ltot, wSuperWeighted, resFull, e_width, resCopy, minFlow ):
	for e in eShort:
		L = e_length[e]
		Ltot += L
		wSuperWeighted += L*resFull[e]/e_width[e]

		resC = resCopy[e]
		if resC < minFlow:
		    minFlow = resC

	wSuperWeighted /= minFlow*Ltot

	return  minFlow, Ltot, wSuperWeighted

def getMinFlowAndWeightsForIntersectPath( eStart, eShort, e_length, Ltot, wSuperWeighted, resFull, e_width, resCopy, minFlow ):
	# AH: There may be a bug here.
	# It is unclear to me why this differs from HSPM version.
	resOld = resFull[eStart]
	for e in eShort:
		L = e_length[e]
		Ltot += L
		wSuperWeighted += L*resOld/e_width[e]
		resOld = resFull[e]            

		resC = resCopy[e]
		if resC < minFlow:
			minFlow = resC

	wSuperWeighted /= minFlow*Ltot

	return  minFlow, Ltot, wSuperWeighted




def getGravitationalAndDiagonalPathParts( vShort, cent, gravPart ):
	v_in   = vShort[1]
	v_out  = vShort[-2]
	zIn    = cent[v_in][2]
	zOut   = cent[v_out][2]
	dz     = zOut -zIn
	pGrav  = 1.0e6 - gravPart*dz

	########
	# change this when using a different DFN size
	dxPart   = (10.6)**2
	########
	Ldiag = np.sqrt( dxPart + (cent[v_out][1] - cent[v_in][1])**2 + dz**2) # radial distance between inflow and outflow

	return pGrav, Ldiag

def getGravitationalPart( vShort, cent, gravPart ):
	# TODO
	# Find out why ISPM requires int here.
	v_in   = int(vShort[1])
	v_out  = int(vShort[-2])
	zIn    = cent[v_in][2]
	zOut   = cent[v_out][2]
	dz     = zOut -zIn
	pGrav  = 1.0e6 - gravPart*dz

	return pGrav


def removeFlowFromPath( eStart, eEnd, eShort, resCopy, minFlow ):
	resCopy[eStart] -= minFlow
	resCopy[eEnd]   -= minFlow
	# remove flow on this path from graph
	for e in eShort:
		resCopy[e] -= minFlow

	return resCopy



def addQforPath( aCmC, wSuperWeighted, Ldiag, pGrav, QgravTot ):
	QgravTot += aCmC/(wSuperWeighted*Ldiag)*pGrav

	return QgravTot