#! /usr/bin/python

# Copyright:
# 	Alex Hobe
# 	Daniel Vogler
# 	Martin P. Seybold


# TODO
# - reduce to functions to see which pieces need to be packed together
# - pack outputs into appropriate classes to reduce passing of arguments to a single object
# - identify when the graph should be initialized for using it effectively as a container


import numpy as np
import copy as cp
import sys
from getFracGeometries import *
from graphMaking import *
from flowRateEstimation import *
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


    graphType = "HananGrid"
    g, aCmC, nVertices, nEdges = makeGraph( intCoordx, intCoordy, intCoordz, graphType, nBoxes, aperture, diff, DX ) 

    cent = g.vertex_properties["cent"]
    e_length      = g.edge_properties["e_length"]
    e_width       = g.edge_properties["e_width"]
    cap           = g.edge_properties["cap"]      # capacity for flow in each edge
    path_crit     = g.edge_properties["path_crit"]

    ############################################################################################
    ## Graph interpretation starts here
    ############################################################################################
    g, max_flow, res, src, tgt = solveMaxFlow(g)

    # terminate if no percolation
    if abs(max_flow - 0.0)<eps:
      print( "\nGraph not connected. No percolation possible.\n")
      return 0.0, 0, nVertices, nEdges # Qvalue, number of paths.

    removeNoFlowEdges(g, res, eps)

    g, res = singleDirectionFlowCleanup( g, res )
    
    removeNoFlowEdges(g, res, eps)

    ############################################################################################
    ### Get path information
    ############################################################################################

    QgravTot, numPaths   = pathFinder( g, src, tgt, res, path_crit, aCmC, eps )
    
    # solver = "HSPM"
    # QgravTot, numPaths, nVertices, nEdges = estimateFlowRate( g, solver )

    return QgravTot, numPaths, nVertices, nEdges
