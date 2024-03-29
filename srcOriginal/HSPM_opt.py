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

def segmentation(workingDir, inputFile):

    tRead = time.time()

    xmldoc     = minidom.parse( workingDir + inputFile )
    # Minidom is used here to cut the xml file into pieces (parsing), add/change data, and creating a new .gz file.
    eps        = 1.0e-13 # used to check for machine precision errors.
    aperture   = 1.0e-5
    diff       = 0.01 - aperture  # used to adjust the arbitrary overlap to capture the nodes in GEOS.

    # retrieve list of matching tags
    itemlist = xmldoc.getElementsByTagName('Nodeset')
    # Inside of the xml file the nodesets are given as follows:
    # <Nodeset name="inlet" type="0"
    #          xmin="-0.01 -0.01 1.99"
    #          xmax="0.01 2.01 2.01" >
    # With these 3 pairs of coordinates all 8 cornerpoints of a fracture are described.
    # The code above thus extracts these fracture coordinates from the xml file.

    # These commands initialize lists and the nBoxes counter which will be filled with the relevant information.
    nBoxes   = 0
    # boxes    = []
    # boxNames = []

    intCoordx = []  # holds the minimum and maximum x-coordinates of the segments
    intCoordy = []  # holds the minimum and maximum y-coordinates of the segments
    intCoordz = []  # holds the minimum and maximum z-coordinates of the segments

    for s in itemlist:
    # Here we run through each fracture (nodeset) in itemlist.
        x,y,z = s.attributes['xmin'].value.split()      # The minimum coordinates are extracted and
        intCoordx.append(float(x) + diff)                             # added to the coordinate lists
        intCoordy.append(float(y) + diff)
        intCoordz.append(float(z) + diff)
        x,y,z = s.attributes['xmax'].value.split()      # same for maximum coordinates
        intCoordx.append(float(x) - diff)
        intCoordy.append(float(y) - diff)
        intCoordz.append(float(z) - diff)
        nBoxes  += 1
    diff = 2*aperture # needed after this point
    # diff = 0.02



    readMeshList = 0
    if readMeshList == 1:
        # Read in domain size and discretization.
        meshList = xmldoc.getElementsByTagName('Mesh')
        for m in meshList:
            xMinDomain, xMaxDomain = m.attributes['xcoords'].value.split()
            yMinDomain, yMaxDomain = m.attributes['ycoords'].value.split()
            zMinDomain, zMaxDomain = m.attributes['zcoords'].value.split()
            nx                     = m.attributes['nx'].value
            ny                     = m.attributes['ny'].value
            nz                     = m.attributes['nz'].value
            DX                     = (float(xMaxDomain) - float(xMinDomain))/int(nx)
            DY                     = (float(yMaxDomain) - float(yMinDomain))/int(ny)
            DZ                     = (float(zMaxDomain) - float(zMinDomain))/int(nz)

    DX = 0.2
    tRead = time.time()- tRead

    tGraph = time.time()
    ############################################################################################
    ### Segmenting the entire domain using the outer coordinates of the fracture intersections
    ### Many of these segments will not correspond to a fracture.
    ############################################################################################

    tG1 = time.time()

    # Reduce intersect coordinate lists to unique coordinate values
    xCoord = np.unique(intCoordx)
    yCoord = np.unique(intCoordy)
    zCoord = np.unique(intCoordz)


    # Initiate interval ID's
    xID = []
    yID = []
    zID = []
    for i in xrange(len(xCoord)-1):
        xID.append([])

    for i in xrange(len(yCoord)-1):
        yID.append([])

    for i in xrange(len(zCoord)-1):
        zID.append([])

    # Find out which fractures lie on which intervals.
    for i in xrange(2,nBoxes):                                # Go through all boxes

        # X-coordinate intervals
        indLow = np.where(xCoord == intCoordx[2*i])         # Find the index of the first interval using xMin
        indUp  = np.where(xCoord == intCoordx[2*i+1])       # Find the index of the last  interval using xMax
        for j in xrange(indLow[0], indUp[0]):     # Add the current fracture number to all intervals between xMin and xMax
            xID[j].append(i)

        # Y-coordinate intervals
        indLow = np.where(yCoord == intCoordy[2*i])         # Find the index of the first interval using yMin
        indUp  = np.where(yCoord == intCoordy[2*i+1])       # Find the index of the last  interval using yMax
        for j in xrange(indLow[0], indUp[0]):     # Add the current fracture number to all intervals between yMin and yMax
            yID[j].append(i)

        # Z-coordinate intervals
        indLow = np.where(zCoord == intCoordz[2*i])         # Find the index of the first interval using yMin
        indUp  = np.where(zCoord == intCoordz[2*i+1])       # Find the index of the last  interval using yMax
        for j in xrange(indLow[0], indUp[0]):     # Add the current fracture number to all intervals between yMin and yMax
            zID[j].append(i)

    tG1 = time.time() - tG1
    tG2 = time.time()
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

    tG2 = time.time() - tG2
    tG3 = time.time()

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

    for k in xrange(0,len(zCoord)-1):
        for j in xrange(0,len(yCoord)-1):
            for i in xrange(0,len(xCoord)-1):
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

    tG3 = time.time() - tG3
    tG4 = time.time()

    ############################################################################################
    ## Attach inlet and outlet to their respective intersections.
    ############################################################################################
    # only one direction is needed
    for i in xrange(0,2):
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

    tG4 = time.time() - tG4

    tGraph = time.time() -tGraph
    ############################################################################################
    ## Graph interpretation starts here
    ############################################################################################

    tMaxFlow = time.time()

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

    tCurrent = time.time()

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

    tPath = time.time() - tPath

    return QgravTot, numPaths, nVertices, nEdges, tRead, tG1, tG2, tG3, tG4, tGraph, tMaxFlow, tBugFix, tPath
