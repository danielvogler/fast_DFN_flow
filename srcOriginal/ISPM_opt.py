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

from xml.dom import minidom
import box3D
import numpy as np
import graph_tool.all as gt # https://graph-tool.skewed.de/static/doc/index.html
import graph_tool.topology as tp
import copy as cp
import sys
import time

#@profile
def simpleMethod(workingDir, inputFile):
    tRead = time.time()    
    xmldoc     = minidom.parse( workingDir + inputFile ) 
    # Minidom is used here to cut the xml file into pieces (parsing), add/change data, and creating a new .gz file.
    
    # retrieve list of matching tags
    itemlist = xmldoc.getElementsByTagName('Nodeset')
    # Inside of the xml file the nodesets are given as follows:
    # <Nodeset name="inlet" type="0"     
    #          xmin="-0.01 -0.01 1.99"     
    #          xmax="0.01 2.01 2.01" >
    # With these 3 pairs of coordinates all 8 cornerpoints of a fracture are described. 
    # The code above thus extracts these fracture coordinates from the xml file.
    
    eps = 1.0e-9 # used for float comparison.
    aperture   = 1.0e-5
    
    # These commands initialize lists and the nBoxes counter which will be filled with the relevant information.
    nBoxes   = 0
    xCoord  = [[],[]] # fracture coordinates
    yCoord  = [[],[]]
    zCoord  = [[],[]]
    verticesInBox = [] # lists which intersects lie on a fracture.
    widthArray = []
    longDimArray = []
    dxArray = []
    dyArray = []
    dzArray = []
    
    for s in itemlist:
    # Here we run through each fracture (nodeset) in itemlist.
        x,y,z = s.attributes['xmin'].value.split()      # The minimum coordinates are extracted and 
        xMin = float(x)
        yMin = float(y)
        zMin = float(z)
        xCoord[0].append(xMin)
        yCoord[0].append(yMin)
        zCoord[0].append(zMin)
        x,y,z = s.attributes['xmax'].value.split()      # same for maximum coordinates
        xMax = float(x)
        yMax = float(y)
        zMax = float(z)
        xCoord[1].append(xMax)
        yCoord[1].append(yMax)
        zCoord[1].append(zMax)
        dx = xMax - xMin
        dy = yMax - yMin
        dz = zMax - zMin
        vals = [dx,dy,dz]
        index = np.argsort(vals) # sorts from smallest to largest and returns index
        widthArray.append(vals[index[1]]) # second largest is the width
        if   dx >= dy and dx >= dz:
            lDim = 0
        elif dy >= dx and dy >= dz:
            lDim = 1
        else:
            lDim = 2
        longDimArray.append(lDim)
        verticesInBox.append( [] )    #init 
        nBoxes  += 1
            
    tRead = time.time()- tRead
    tGraph = time.time()
    tG1 = time.time()
        
    # intersect boxes to retain vertices 
    vertexPos  = []
    vertexPos.append([-0.3,5,5])
    vertexPos.append([10.3,5,5])
 

    # build graph
    g = gt.Graph()  # Initialize a new gt graph.
    g.set_directed(True) # directed graphs have separate amounts for flow in each direction.
    

    nVertices  = 2
    
    # calculate edges, edge lengths and box widths
    e_length  = g.new_edge_property("double")
    e_width   = g.new_edge_property("double")
    cap       = g.new_edge_property("double")      # capacity for flow in each edge
    path_crit = g.new_edge_property("double")      # criterion for finding the shortest path
    
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
                vertexPos.append( [xCent, yCent, zCent] ) # saves the ceter point (node) coordinates
                nVertices +=1
    g.add_vertex(nVertices)   # Initializing all intersection vertices
    
    tG1 = time.time() - tG1
    tG2 = time.time()     
        

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
 
    tG2 = time.time() - tG2
    tG3 = time.time()
    
    # Normalize to avoid numerical precision issues.
    capIndex = np.where(cap.a < 99)
    cap_max = max(cap.a[capIndex])
    cap.a[capIndex] /= cap_max


    tG3 = time.time() - tG3     
    tG4 = time.time()
    

    ############################################################################################
    ## Attach inlet and outlet to their respective intersections.
    ############################################################################################
    # only one direction is needed
    for i in xrange(0,2):
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

    tG4 = time.time() - tG4

    # print "Graph time: ", time.time() -tGraph
    tGraph = time.time() -tGraph
    tMaxFlow = time.time() 
    
    # solving the max flow problem
    full_source = 0
    full_target = 1
    src, tgt = g.vertex(full_source), g.vertex(full_target)               # set the source and target nodes
    res = gt.edmonds_karp_max_flow(g, src, tgt, cap)  
    res.a = cap.a - res.a  # the actual flow
    max_flow = sum(res[e] for e in tgt.in_edges())    # the maximum flow is the sum of the individual paths

    if abs(max_flow)<eps:
      	print ""
      	print "Graph not connected. No percolation possible."
      	print ""
      	return 0.0, 0, nVertices, nEdges # Qvalue, number of paths.
    
    noFlow = gt.find_edge_range(g, res, [0.0, eps])
    for s in noFlow:
      g.remove_edge(s)

    tMaxFlow = time.time() - tMaxFlow

    
    tBugFix = time.time()
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

    tBugFix = time.time() -tBugFix
    
    
    
    
    tPath = time.time()
    
    ############################################################################################
    ### Get path information
    ############################################################################################
    
    def pathFinder(gCopy, src, tgt, resCopy, path_crit, aCmC):

      numPaths     = 0
      resFull      = resCopy.copy()
      QgravTot     = 0.
      gravPart = 9.81*1000.
      
      for i in xrange(0,1000): # perhaps do a while loop instead. This should be safer...
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
      
          wSuperWeighted = 0.0
          resOld = resFull[eStart]
          for e in eShort:
                L = e_length[e]
                wSuperWeighted += L*resOld/e_width[e]
                resOld = resFull[e]            
                resC = resCopy[e]
                if resC < minFlow:
                  minFlow = resC
          wSuperWeighted /= minFlow

          v_in   = int(vShort[1])
          v_out  = int(vShort[-2])
          pGrav  = 1.0e6 - gravPart*(vertexPos[v_out][2] - vertexPos[v_in][2])

          resCopy[eStart] -= minFlow
          resCopy[eEnd]   -= minFlow
          for i,e in enumerate(eShort):
            # remove flow on this path from graph
            resCopy[e] -= minFlow
            
          # calculate Q for this path
          QgravTot += aCmC/(wSuperWeighted)*pGrav
          
          # remove edges without flow 
          noFlow = gt.find_edge_range(gCopy, resCopy, [0.0, eps])
          for s in noFlow:
            gCopy.remove_edge(s)
    
          numPaths += 1
    
    QgravTot, numPaths = pathFinder(g, src, tgt, res, path_crit, aCmC)

    tPath = time.time() - tPath

    return QgravTot, numPaths, nVertices, nEdges, tRead, tG1, tG2, tG3, tG4, tGraph, tMaxFlow, tBugFix, tPath
    

