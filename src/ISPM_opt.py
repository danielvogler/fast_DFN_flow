#! /usr/bin/python

# Copyright:
# 	Alex Hobe
# 	Daniel Vogler
# 	Martin P. Seybold

from xml.dom import minidom
import numpy as np
# import graph_tool.all as gt # https://graph-tool.skewed.de/static/doc/index.html
# import graph_tool.topology as tp
import copy as cp
import sys
from getFracGeometries import *
from graphMaking import *
from flowRateEstimation import *
from graphToolwrapper import *

gt = graphToolwrapper.all()
tp = graphToolwrapper.topology()

#@profile
def simpleMethod(workingDir, inputFile):
    eps = 1.0e-9 # used for float comparison.
    aperture   = 1.0e-5

    itemlist = getCoordListGEOSxml(workingDir + inputFile)   

    xCoord, yCoord, zCoord, nBoxes, longDimArray, verticesInBox, widthArray = makeFracBoxGEOSxmlForISPM( itemlist )
            
    # intersect boxes to retain vertices 
    vertexPos  = []
    vertexPos.append([-0.3,5,5])
    vertexPos.append([10.3,5,5])
 

    # build graph
    # g = gt.Graph()  # Initialize a new gt graph.
    # g.set_directed(True) # directed graphs have separate amounts for flow in each direction.
    sourceTargetCentroids = [   [-0.3, 5, 5],
                                [10.3, 5, 5]
                            ]
    g = initializeGraph( sourceTargetCentroids )
    

    nVertices  = 2

    properties = [  ["e_length", "e_width", "cap", "path_crit"],
                    ["double", "double", "double", "double"]
                  ]
    g = makeNewEdgeProperty( g, properties )
    e_length      = g.edge_properties["e_length"]
    e_width       = g.edge_properties["e_width"]
    cap           = g.edge_properties["cap"]      # capacity for flow in each edge
    path_crit     = g.edge_properties["path_crit"]
    
    # calculate edges, edge lengths and box widths
    
    for i in range( 0, nBoxes ):
        xFracMin = xCoord[0][i]
        yFracMin = yCoord[0][i]
        zFracMin = zCoord[0][i]
        xFracMax = xCoord[1][i]
        yFracMax = yCoord[1][i]
        zFracMax = zCoord[1][i]
        for j in range(i+1, nBoxes):
            xIntMin = max(xFracMin,xCoord[0][j])
            xIntMax = min(xFracMax,xCoord[1][j])
            yIntMin = max(yFracMin,yCoord[0][j])
            yIntMax = min(yFracMax,yCoord[1][j])
            zIntMin = max(zFracMin,zCoord[0][j])
            zIntMax = min(zFracMax,zCoord[1][j])
            if (xIntMax > xIntMin) and (yIntMax > yIntMin) and (zIntMax > zIntMin):
                verticesInBox[i].append(nVertices) # saves the two
                verticesInBox[j].append(nVertices) # fractures that intersect
                xCent = 0.5*(xIntMax+xIntMin)
                yCent = 0.5*(yIntMax+yIntMin)
                zCent = 0.5*(zIntMax+zIntMin)
                vertexPos.append( [xCent, yCent, zCent] ) # saves the center point (node) coordinates
                nVertices +=1
    g.add_vertex(nVertices)   # Initializing all intersection vertices
    

    nEdges = 0 # count edge number    
    mu            = 0.001
    cubicLawConst = 12
    aCmC          = aperture**3/(mu*cubicLawConst)
    cap_max = 0
    vertexPos = np.array(vertexPos)
    eprops = [e_width, e_length, path_crit, cap]
    for j in range(2, nBoxes): # Go through all fractures while leaving out the inlet and outlet boxes. These will be added later. 
        vInBoxSorted = verticesInBox[j] 
        nVinFrac = len(vInBoxSorted) -1
        if nVinFrac <1:
            continue
        lDim  = longDimArray[j]
        vInBoxSorted = sorted( vInBoxSorted, key = lambda x : vertexPos[x][lDim] )
        widthList = (widthArray[j]-0.02)*np.ones(2*nVinFrac)
        widthList2 = (widthArray[j]-0.02)*np.ones(nVinFrac)
        edgeList = [[],[],[],[],[],[]]
        v1 = vInBoxSorted[:-1]
        v2 = vInBoxSorted[1:]
        Vv1 = v1+v2
        Vv1[::2] = v1
        Vv1[1::2] = v2
        Vv2 = v2+v1
        Vv2[::2]  = v2
        Vv2[1::2] = v1
        edgeList[0] = Vv1
        edgeList[1] = Vv2
        edgeList[2] = widthList
        p1x = vertexPos[v1][:,0]
        p1y = vertexPos[v1][:,1]
        p1z = vertexPos[v1][:,2]
        p2x = vertexPos[v2][:,0]
        p2y = vertexPos[v2][:,1]
        p2z = vertexPos[v2][:,2]
        lengthArray = ( (p1x - p2x)**2 + (p1y - p2y)**2 + (p1z - p2z)**2 )**0.5
        doubleLengthArray = np.repeat(lengthArray,2)
        edgeList[3] = doubleLengthArray
        instantFlow = np.where(lengthArray <= 1.e-9)
        edgeList[4] = doubleLengthArray/widthList

        lengthArray[instantFlow] = 1.e-9
        capArray =  aCmC*widthList2/lengthArray
        capArray[instantFlow] = 100.
        edgeList[5] = np.repeat(capArray,2)
        g.add_edge_list(zip(*edgeList), eprops=eprops)
        nEdges += 2*nVinFrac
 
    # Normalize to avoid numerical precision issues.
    cap, cap_max = normalizeEdgeCapacities(g)

    ############################################################################################
    ## Attach inlet and outlet to their respective intersections.
    ############################################################################################
    # only one direction is needed
    for i in range(0,2):
        vInBox = verticesInBox[i]
        lV     = len(vInBox)
        oneArray = np.ones(lV)
        bigV   = i*oneArray
        lengthArray = 0*oneArray
        capArray  = 100*oneArray
        edgeList = [[],[],[],[],[],[]]

        if i == 0:
            edgeList[0] = bigV
            edgeList[1] = vInBox
        else:
            edgeList[0] = vInBox
            edgeList[1] = bigV

        edgeList[2] = lengthArray # widths are zero
        edgeList[3] = lengthArray # lengths are zero
        edgeList[4] = lengthArray # path_crits are zero
        edgeList[5] = capArray 
        g.add_edge_list(zip(*edgeList), eprops=eprops)
        nEdges += lV

    ############################################################################################
    ## Graph interpretation starts here
    ############################################################################################ 
    
    # solving the max flow problem
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
    
    def pathFinder(gCopy, src, tgt, resCopy, path_crit, aCmC):

      numPaths, resFull, QgravTot, gravPart = initializePathFinder(resCopy)
      
      for i in range(0,1000): # perhaps do a while loop instead. This should be safer...

        
        vShort, eShort = tp.shortest_path(gCopy,src,tgt, weights=path_crit) # find shortest path
        
        if len(eShort) == 0: # no more path found
          return QgravTot, numPaths  # return values
        
        else:
          pLength, pWidth, pFactor, pRes, minFlow, wSuperWeighted, Ltot = initializeIndividualPath()

          eStart, eEnd, eShort =  getAndRemoveStartAndEndFromPath(eShort)

          # minFlow, Ltot, wSuperWeighted = getMinFlowAndWeightsForHananPath( eShort, e_length, Ltot, wSuperWeighted, resFull, e_width, resCopy, minFlow )

 
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

          pGrav = getGravitationalPart( vShort, vertexPos, gravPart )

          resCopy = removeFlowFromPath( eStart, eEnd, eShort, resCopy, minFlow )
            
          QgravTot = addQforPath( aCmC, wSuperWeighted, Ltot, pGrav, QgravTot )

          
          removeNoFlowEdges(gCopy, resCopy, eps)
    
          numPaths += 1
    
    QgravTot, numPaths = pathFinder(g, src, tgt, res, path_crit, aCmC)
    # QgravTot, numPaths   = pathFinder( g, src, tgt, res, path_crit, aCmC, eps )

    return QgravTot, numPaths, nVertices, nEdges
    

