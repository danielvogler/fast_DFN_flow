#! /usr/bin/python

# Copyright:
# 	Alex Hobe
# 	Daniel Vogler
# 	Martin P. Seybold


# TODO
# - pack outputs into appropriate classes to reduce passing to one object



from xml.dom import minidom
import numpy as np
import graph_tool.all as gt # https://graph-tool.skewed.de/static/doc/index.html
import graph_tool.topology as tp
import copy as cp
import sys
from getFracGeometries import *

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

    # Reduce intersect coordinate lists to unique coordinate values
    xCoord = np.unique(intCoordx)
    yCoord = np.unique(intCoordy)
    zCoord = np.unique(intCoordz)


    # Initiate interval ID's
    xID = []
    yID = []
    zID = []
    for i in range(len(xCoord)-1):
        xID.append([])

    for i in range(len(yCoord)-1):
        yID.append([])

    for i in range(len(zCoord)-1):
        zID.append([])

    # Find out which fractures lie on which intervals.
    for i in range(2,nBoxes):                                # Go through all boxes

        # X-coordinate intervals
        indLow = np.where(xCoord == intCoordx[2*i])         # Find the index of the first interval using xMin
        indUp  = np.where(xCoord == intCoordx[2*i+1])       # Find the index of the last  interval using xMax
        for j in range(indLow[0][0], indUp[0][0]):     # Add the current fracture number to all intervals between xMin and xMax
            xID[j].append(i)

        # Y-coordinate intervals
        indLow = np.where(yCoord == intCoordy[2*i])         # Find the index of the first interval using yMin
        indUp  = np.where(yCoord == intCoordy[2*i+1])       # Find the index of the last  interval using yMax
        for j in range(indLow[0][0], indUp[0][0]):     # Add the current fracture number to all intervals between yMin and yMax
            yID[j].append(i)

        # Z-coordinate intervals
        indLow = np.where(zCoord == intCoordz[2*i])         # Find the index of the first interval using yMin
        indUp  = np.where(zCoord == intCoordz[2*i+1])       # Find the index of the last  interval using yMax
        for j in range(indLow[0][0], indUp[0][0]):     # Add the current fracture number to all intervals between yMin and yMax
            zID[j].append(i)

    ############################################################################################
    ### graph creation
    ############################################################################################
    g = gt.Graph()          # Initialize a new gt graph.
    g.set_directed(True)    # directed graphs have separate amounts for flow in each direction.

    cent  = g.new_vertex_property("vector<double>")   # node (centroid) coordinates for path positioning.

    # define source of the graph
    vI = g.add_vertex()
    cent[ vI ]  = [0.0,0.0,0.0]

    # define target of the graph
    vI = g.add_vertex()
    cent[ vI ]  = [0.0,0.0,0.0]

    # find the segments that lie on fractures and turn them into nodes.
    nVertices  = 0
    domainSegBoxes = [[],[],[]]
    segments = -1*np.ones((len(zCoord)-1, len(yCoord)-1,len(xCoord)-1))

    vertsInBox = [[],[]] # lists which segments lie on a fracture

    for i, xIDs in enumerate(xID):                 # sweep through x-coordinate intervals
        for j, yIDs in enumerate(yID):             # find all segments that have lie on this x-coordinate slice
            yResult = set(xIDs) & set(yIDs)
            if not yResult:                        # an empty set is False.
                continue
            else:
                for k,zIDs in enumerate(zID):
                    zResult = yResult & set(zIDs)
                    if not zResult:                        # an empty set is False.
                        continue
                    else:
                        xMin  = xCoord[i]
                        xMax  = xCoord[i+1]
                        yMin  = yCoord[j]
                        yMax  = yCoord[j+1]
                        zMin  = zCoord[k]
                        zMax  = zCoord[k+1]
                        domainSegBoxes[0].append(xMax - xMin)
                        domainSegBoxes[1].append(yMax - yMin)
                        domainSegBoxes[2].append(zMax - zMin)

                        xCent = 0.5*(xMax + xMin)
                        yCent = 0.5*(yMax + yMin)
                        zCent = 0.5*(zMax + zMin)

                        segments[k][j][i] = nVertices
                        vI = g.add_vertex()
                        cent[ vI ]  = [xCent, yCent, zCent]

                        if i == 0:
                            vertsInBox[0].append(nVertices)
                        elif i == len(xID)-1:
                            vertsInBox[1].append(nVertices)
                        nVertices +=1
    nVertices += 2


    ############################################################################################
    ### Connecting domain segments that lie inside a fracture
    ############################################################################################

    e_length      = g.new_edge_property("double")
    e_width       = g.new_edge_property("double")
    cap           = g.new_edge_property("double")      # capacity for flow in each edge
    path_crit     = g.new_edge_property("double")      # criterion for finding the shortest path
    nEdges        = 0
    mu            = 0.001
    cubicLawConst = 12
    aCmC          = aperture**3/(mu*cubicLawConst)

    def connect_nodes(node1,node2,width, diff, DX, nEdges, aCmC):
        E_ps = 1.0e-9
        # This part ignores the arbitrary value in the xml file by not connecting edges with a width of diff m
        # and by extending or reducing those that are slightly too small or slightly too large.
        if abs(width - diff) <1.0e-13:       # Find edges with width == diff m
            return nEdges # ignore this edge


        v1, v2 = g.vertex(node1+2) ,  g.vertex(node2+2)         # Get the vertex labels for this vertex and the next vertex for this fracture.
        p1, p2 = cent[v1],cent[v2]   # Get the position of these two vertices.
        e = g.add_edge(v1, v2)                                  # Create an edge between these two vertices.
        eL = ( abs(p1[0] - p2[0])**2 + abs(p1[1] - p2[1])**2 + abs(p1[2] - p2[2])**2 )**.5 # Length of the edge calculated as the absolute distance between two vertices.
        e_length[e] = eL
        e_width[e]  = width                                     # this width corresponds to the width direction
        assert e_length[e] >= 0.0 # This checks if the assumption of a positive edge length is correct.

        #These following commands create the edge in the opposite direction.
        eBack = g.add_edge(v2, v1)
        e_length[eBack] = eL
        e_width[eBack]  = width

        if eL <= 1e-9:                   # set capacity for overlaping nodes.
            cap[e]           = 100.0
            cap[eBack]       = 100.0
            path_crit[e]     =   0.0                     # criterion for finding the shortest path
            path_crit[eBack] =   0.0                     # criterion for finding the shortest path

        else :
            eLeW             = eL/width
            eCap             = aCmC / eLeW
            cap[e]           = eCap
            cap[eBack]       = eCap
            path_crit[e]     = eLeW                 # criterion for finding the shortest path
            path_crit[eBack] = eLeW                 # criterion for finding the shortest path

        nEdges += 2
        return nEdges

    ####
    # Go through the domain and attach the nodes with edges.
    ####

    for k in range(0,len(zCoord)-1):
        for j in range(0,len(yCoord)-1):
            for i in range(0,len(xCoord)-1):
                seg = int(segments[k][j][i])  # get segment
                if seg== -1:
                    continue
                else:
                    dx = domainSegBoxes[0][seg]
                    dy = domainSegBoxes[1][seg]
                    dz = domainSegBoxes[2][seg]

                    # check neighbours
                    if i<len(xCoord)-2:
                        segright = int(segments[k][j][i+1])
                        if segright !=-1: # connect these
                            width = max(dy,dz)
                            nEdges = connect_nodes(seg,segright,width, diff, DX, nEdges, aCmC)

                    if j<len(yCoord)-2:
                        segback  = int(segments[k][j+1][i])
                        if segback !=-1: # connect these
                            width = max(dx,dz)
                            nEdges = connect_nodes(seg,segback,width, diff, DX, nEdges, aCmC)

                    if k<len(zCoord)-2:
                        segup    = int(segments[k+1][j][i])
                        if segup !=-1: # connect these
                            width = max(dx,dy)
                            nEdges = connect_nodes(seg,segup,width, diff, DX, nEdges, aCmC)

    # Normalize to avoid numerical precision issues.
    capIndex = np.where(cap.a < 99)
    cap_max = max(cap.a[capIndex])
    cap.a[capIndex] /= cap_max


    ############################################################################################
    ## Attach inlet and outlet to their respective intersections.
    ############################################################################################
    # only one direction is needed
    for i in range(0,2):
        for j in vertsInBox[i]:
            if i == 0:
                v1, v2 = g.vertex(i) ,  g.vertex(j+2) # source to boundary
            else:
                v1, v2 = g.vertex(j+2), g.vertex(i)   # boundary to target

            e = g.add_edge(v1, v2)

            e_length[e] = 0.0
            e_width[e]  = 0.0
            assert e_length[e] >= 0.0
            cap[e]       = 100.0
            path_crit[e] = 0.0                     # criterion for finding the shortest path
            nEdges += 1

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
