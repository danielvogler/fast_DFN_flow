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

	if solver == "SPM": # shortest path maxflow
		g, max_flow = solveMaxFlow(cap)

		# terminate of no percolation
		if abs(max_flow - 0.0)<eps:
			print( "" )
			print( "Graph not connected. No percolation possible." )
			print( "" )
			return 0.0, 0, nVertices, nEdges # Qvalue, number of paths.

		removeNoFlowEdges(g, res, eps)
		# Cleanup before next step to speed things up

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

		removeNoFlowEdges(g, res, eps)
		# remaining Cleanup

		QgravTot, numPaths   = pathFinder(g, src, tgt, res, path_crit, aCmC)
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





def pathFinder(gCopy, src, tgt, resCopy, path_crit, aCmC):
	# TODO
	# - fix indentation
	# - test resCopy.copy()


	numPaths     = 0
	resFull      = resCopy.copy()
	QgravTot     = 0.

	gravPart = 9.81*1000
	########
	# change this when using a different DFN size
	dxPart   = (10.6)**2
	########

	for i in xrange(0,1000): # perhaps do a while loop instead. This should be safer...
		pLength    = []  # store lengths of this path
		pWidth     = []  # store widths of this path
		pFactor    = []  # store flow factor of this path (percentage of the edge that this path uses).
		pRes       = []  # store the flow results for this path


		vShort, eShort = tp.shortest_path(gCopy,src,tgt, weights=path_crit) # find shortest path

		if len(eShort) == 0: # no more path found
		    return QgravTot, numPaths

		else:
			minFlow = 10000

			eStart = eShort[0]
			eEnd   = eShort[-1]
			del eShort[0]
			del eShort[-1]

			Ltot = 0.
			wSuperWeighted = 0.0
			for e in eShort:
				L = e_length[e]
				Ltot += L
				wSuperWeighted += L*resFull[e]/e_width[e]

				resC = resCopy[e]
				if resC < minFlow:
				    minFlow = resC

			wSuperWeighted /= minFlow*Ltot


			v_in   = vShort[1]
			v_out  = vShort[-2]
			zIn    = cent[v_in][2]
			zOut   = cent[v_out][2]
			dz     = zOut -zIn
			pLdiag = np.sqrt( dxPart + (cent[v_out][1] - cent[v_in][1])**2 + dz**2) # radial distance between inflow and outflow
			pGrav  = 1.0e6 - gravPart*dz


			resCopy[eStart] -= minFlow
			resCopy[eEnd]   -= minFlow
			# remove flow on this path from graph
			for e in eShort:
				resCopy[e] -= minFlow

			# calculate Q for this path
			QgravTot += aCmC/(wSuperWeighted*pLdiag)*pGrav

			removeNoFlowEdges(g, res, eps)
			numPaths += 1
