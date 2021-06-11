#! /usr/bin/python

# Copyright:
# 	Alex Hobe
# 	Daniel Vogler
# 	Martin P. Seybold


# TODO
# - adjust to using fracBoxes
# - abstract should take fracGeoms as input instead of fracBoxes

import sys # import sys.exit only?
import numpy as np
# import graph_tool.all as gt # https://graph-tool.skewed.de/static/doc/index.html

from graphToolwrapper import *

gt = graphToolwrapper.all()



def abstractGraphMaker(fracBoxes, graphType):
	""" Set of methods, that create a graph out of the fracture geometries.
	"""

def makeGraph(abstractGraphMaker):
	# TODO
	# -rename graphType
	# 	- nodeEdgeConvention
	# 	- nodeEdgeAssumption
	# - extract common parts from if clause
	# - use inheritance instead of hardcoded if clause.
	# - rename cent > nodeCentroidCoords?
	# - how to pass around diff, DX?
	# - reduce output to one object
	# 	- are new_edge_properties assessible using g.property?
	#

	if graphType == "HananGrid":
		xID, yID, zID = segmentDomain( intCoordx, intCoordy, intCoordz )
		g = initializeGraph( )
		cent, vertsInBox, segments, domainSegBoxes = addFracSegmentNode( g,xID, yID, zID,xCoord,yCoord,zCoord )
		e_length, e_width, cap, path_crit, nEdges  = addEdgesBetweenFracSegments( g,xCoord, yCoord, zCoord, segments, domainSegBoxes )
		
	elif graphType == "IntersectSize":
		print( " graphType: 'IntersectSize' not yet implemented" )
		# TODO: error/exception instead of print
		sys.exit()

	
	# return fracGraph
	return g

def segmentDomain( intCoordx, intCoordy, intCoordz ):
	# TODO
	# - think about abstraction of this after refactoring ISPM version.
	# - refactor further after testing
	
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
	
	return xID, yID, zID


def initializeGraph( ):
	# TODO
	# - change centroid position of target?

	g = gt.Graph() 
	g.set_directed(True) 
	cent  = g.new_vertex_property("vector<double>")   # node (centroid) coordinates for path positioning.
    
	# define source of the graph
    vI = g.add_vertex()
    cent[ vI ]  = [0.0,0.0,0.0]
    
    # define target of the graph
    vI = g.add_vertex()                                    
    cent[ vI ]  = [0.0,0.0,0.0]

	# Notes:
	# - Directed graphs have separate amounts for flow in each direction.
	# - Do not use set_fast_edge_removal! This changes edge ordering.
	# 
	
	return g

def addFracSegmentNode( g,xID, yID, zID, xCoord, yCoord, zCoord ):
	# TODO
	# - think about abstraction of this after refactoring ISPM version.
	# - refactor further after testing

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
    nVertices += 2 # for source and target, which are not in vertsInBox
	
	return cent, vertsInBox, segments, domainSegBoxes

def addEdgesBetweenFracSegments( g,xCoord, yCoord, zCoord, segments, domainSegBoxes ):
	# TODO
	# - think about abstraction of this after refactoring ISPM version.
	# - refactor further after testing

    # initialize 
    e_length      = g.new_edge_property("double")
    e_width       = g.new_edge_property("double")
    cap           = g.new_edge_property("double")      # capacity for flow in each edge
    path_crit     = g.new_edge_property("double")      # criterion for finding the shortest path
    nEdges        = 0
    mu            = 0.001
    cubicLawConst = 12
    aCmC          = aperture**3/(mu*cubicLawConst)

    # Go through the domain and attach the nodes with edges.
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

    # Attach inlet and outlet to their respective intersections.
    e_length, e_width, cap, path_crit, nEdges = attachBoundaryNodes( vertsInBox )

	return e_length, e_width, cap, path_crit, nEdges 

def connect_nodes( node1,node2,width, diff, DX, nEdges, aCmC ):
	E_ps = 1.0e-9
    # This part ignores the arbitrary value in the xml file by not connecting edges 
	# with a width of diff.        
    if abs(width - diff) <1.0e-13:       # Find edges with width == diff
    	return nEdges # ignore this edge

    v1, v2 = g.vertex(node1+2) ,  g.vertex(node2+2)         # Get the vertex labels for this vertex and the next vertex for this fracture.
    p1, p2 = cent[v1],cent[v2]   # Get the position of these two vertices.
    e = g.add_edge(v1, v2)                                  # Create an edge between these two vertices.
    # list which fracture this edge belonged to.
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

def attachBoundaryNodes( vertsInBox ):
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

    return e_length, e_width, cap, path_crit, nEdges




