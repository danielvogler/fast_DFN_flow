#! /usr/bin/python

# Copyright:
# 	Alex Hobe
# 	Daniel Vogler
# 	Martin P. Seybold


# TODO
# - pack outputs into appropriate classes to reduce passing to one object


import numpy as np
import copy as cp
import sys
from getFracGeometries import *
from graphMaking import *
from graphToolwrapper import *

gt = graphToolwrapper.all()
tp = graphToolwrapper.topology()

def segmentation(workingDir, inputFile):


    eps        = 1.0e-13 # used to check for machine precision errors.
    aperture   = 1.0e-5
    
    fileType = "GEOSxml"

    # get Fracture information
    intCoordx, intCoordy, intCoordz, nBoxes = getFractureBoxes( workingDir + inputFile, fileType )


    diff = 2*aperture # needed after this point
    DX = 0.2

    ############################################################################################
    ### Segmenting the entire domain using the outer coordinates of the fracture intersections
    ### Many of these segments will not correspond to a fracture.
    ############################################################################################

    xID, yID, zID, xCoord, yCoord, zCoord = segmentDomain( intCoordx, intCoordy, intCoordz, nBoxes )

    
    ############################################################################################
    ### graph creation
    ############################################################################################
    sourceTargetCentroids = [[0.0,0.0,0.0],[0.0,0.0,0.0]]
    g = initializeGraph( sourceTargetCentroids )


    vertsInBox, segments, domainSegBoxes, nVertices = addFracSegmentNode( g,xID, yID, zID, xCoord, yCoord, zCoord )
    # g is changed here without needing to be returned by the function.
    # this is not "clean code" best practice...
    cent = g.vertex_properties["cent"]  


    e_length, e_width, cap, path_crit, nEdges, aCmC = addEdgesBetweenFracSegments( g,xCoord, yCoord, zCoord, segments, domainSegBoxes, aperture, diff, DX, vertsInBox )

    ############################################################################################
    ## Graph interpretation starts here
    ############################################################################################

    # solving the max flow problem
    full_source = 0
    full_target = 1
    src, tgt = g.vertex(full_source), g.vertex(full_target)               # set the source and target nodes
    res = gt.edmonds_karp_max_flow(g, src, tgt, cap)
    res.a = cap.a - res.a  # the actual flow
    max_flow = sum(res[e] for e in tgt.in_edges())    # the maximum flow is the sum of the individual paths

    if abs(max_flow - 0.0)<eps:
      print( "\nGraph not connected. No percolation possible.\n")
      return 0.0, 0, nVertices, nEdges # Qvalue, number of paths.

    noFlow = gt.find_edge_range(g, res, [0.0, eps])
    for s in noFlow:
        g.remove_edge(s)

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

    noFlow = gt.find_edge_range(g, res, [0.0, eps])
    for s in noFlow:
        g.remove_edge(s)

    ############################################################################################
    ### Get path information
    ############################################################################################

    def pathFinder(gCopy, src, tgt, resCopy, path_crit, aCmC):

      numPaths     = 0
      resFull      = resCopy.copy()
      QgravTot     = 0.

      gravPart = 9.81*1000
      ########
      # change this when using a different DFN size
      dxPart   = (10.6)**2
      ########

      for i in range(0,1000): # perhaps do a while loop instead. This should be safer...
        pLength    = []  # store lengths of this path
        pWidth     = []  # store widths of this path
        pFactor    = []  # store flow factor of this path (percentage of the edge that this path uses).
        pRes       = []  # store the flow results for this path


        vShort, eShort = tp.shortest_path(gCopy,src,tgt, weights=path_crit) # find shortest path

        if len(eShort) == 0: # no more path found
          return QgravTot, numPaths  # return values

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
          for e in eShort:
            # remove flow on this path from graph
            resCopy[e] -= minFlow

          # calculate Q for this path
          QgravTot += aCmC/(wSuperWeighted*pLdiag)*pGrav

          # remove edges without flow
          noFlow = gt.find_edge_range(gCopy, resCopy, [0.0, eps])
          for s in noFlow:
            gCopy.remove_edge(s)
          numPaths += 1

    QgravTot, numPaths   = pathFinder(g, src, tgt, res, path_crit, aCmC)

    return QgravTot, numPaths, nVertices, nEdges
