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


import sys
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colorbar as colorbar
from numpy import *
from xml.dom import minidom

import box3D
import graph_tool.all as gt # https://graph-tool.skewed.de/static/doc/index.html
from overlap_and_boundary import *
from generateMsub import *
import os

def generateDFN(simulation, workingDir, simName, numberOfFracturesNeeded, savePNG, GEOSjob):
    #############################################################################
    ## user input
    #############################################################################
    ##
    ## Description and setting of fracture network parameters.
    # choose parameters for fracture network
    # Desired amount of fractures
    name = simName + "_{0}".format(simulation)     # filename
    overshoot = int(np.round(0.9*numberOfFracturesNeeded))
    numberOfFractures = numberOfFracturesNeeded + overshoot    # Add an overshoot, because overlapping fractures will be removed.
    distribution = 0 # use a uniform distribution between min and max. if 0 us min by max fractures
    #print( "number of fractures needed.", numberOfFracturesNeeded)
    #print( "number of fractures created.", numberOfFractures)

    # mesh partitioning
    nx = 50
    ny = 50    # The fracture network needs to be discretized on a mesh for GEOS to be able to work.
    nz = 50    # Here the amount of nodes in the x, y, and z-directions are defined.

    xSize = 10.0
    ySize = 10.0
    zSize = 10.0

    # box dimensions
    xcoords = array([0.0, xSize])     # The dimension of the box is defined by setting the beginning and end coordinates of each axis. 
    ycoords = array([0.0, ySize])     # The dimension in a particular direction is the second coordinate minus the first in meters.
    zcoords = array([0.0, zSize]) 

    # fracture size parameters
    fractureSizeMin = 2.4    # in meters
    fractureSizeMax = 3.4
    diff = 0.01            # used to ensure GEOS node is within the fracture box
    extend = 0.2           # makes each fracture longer to ensure proper connections are made

    Pin  = 2.0e6           # inlet pressure
    Pout = 1.0e6           # outlet pressure

    ##
    #################################################################################

    # This code will use these parameters to randomly distribute fractures parallel to the x-y-plane, the x-z-plane, and the y-z-plane. 
    # The fractures are rectangles. The length of these rectangles follow a uniform distribution between the minimum and maximum size defined here. 
    # The orientation is a uniform discrete distribution.


    # Create random sized rectangular fractures with random location and  orientation.
    # calculate fracture edge length of the two rectangle sides
    fractureEdgeLength      = zeros((2,numberOfFractures))

    if distribution == 1:
        fractureEdgeLength[0,:] = (fractureSizeMax - fractureSizeMin)*random.random((numberOfFractures,)) + fractureSizeMin + 2*extend
        fractureEdgeLength[1,:] = (fractureSizeMax - fractureSizeMin)*random.random((numberOfFractures,)) + fractureSizeMin + 2*extend
    else:
        fractureEdgeLength[0,:] = fractureSizeMin #+ 2*extend
        fractureEdgeLength[1,:] = fractureSizeMax #+ extend
    # Here the two edge lengths of the rectangular fractures are created using a uniform distribution between the minimum and maximum lengths set above. 

    # create random fracture locations
    fractureLocation        = zeros((3,numberOfFractures))
    fractureLocation[0,:]   = random.randint(1, nx, size=numberOfFractures) * ( xcoords[1] - xcoords[0] )/nx + xcoords[0]
    fractureLocation[1,:]   = random.randint(1, ny, size=numberOfFractures) * ( ycoords[1] - ycoords[0] )/ny + ycoords[0]
    fractureLocation[2,:]   = random.randint(1, nz, size=numberOfFractures) * ( zcoords[1] - zcoords[0] )/nz + zcoords[0]
    # This one creates a uniform discrete distribution of center points. 
    # For each fracture one mid point is defined. 
    # This is then combined with the fracture size described above and the orientation described below to define the final fracture. 
        
    # fracture orientation
    # 1 - [1 0 0]
    # 2 - [0 1 0]
    # 3 - [0 0 1]
    fractureOrientation     = random.randint(1, 4, size=numberOfFractures) # chooses from min = 1 and max = 4-1 (python syntax)
    # a uniform discrete distribution. 
        

    ## Creating the discrete fracture network with the assigned random numbers.
    fractureCornerPoints = zeros((numberOfFractures,4,3))
    for i in range(0,numberOfFractures):
    # This loop finds the four corner points for each fracture according to its orientation and the two side lengths. 
    # The four points are calculated by permutating one coordinate +/- half the edge length in that direction and the other coordinate +/- the edge length in this other direction. 
    # It is thus possible for the fractures to cross the boundaries.

        # normal vector x-direction
        if fractureOrientation[i] == 1:
            edgeLengthFactor = array([0, 1 ,1]) # Does not seem to be used anywhere.
            fractureCornerPoints[i,0,:] = array([ (fractureLocation[0,i] ) , \
                                        (fractureLocation[1,i] + fractureEdgeLength[0,i]/2) ,\
                                        (fractureLocation[2,i] + fractureEdgeLength[1,i]/2) ])
            fractureCornerPoints[i,1,:] = array([ (fractureLocation[0,i] ) , \
                                        (fractureLocation[1,i] + fractureEdgeLength[0,i]/2) ,\
                                        (fractureLocation[2,i] - fractureEdgeLength[1,i]/2) ])
            fractureCornerPoints[i,2,:] = array([ (fractureLocation[0,i] ) , \
                                        (fractureLocation[1,i] - fractureEdgeLength[0,i]/2) ,\
                                        (fractureLocation[2,i] - fractureEdgeLength[1,i]/2) ])
            fractureCornerPoints[i,3,:] = array([ (fractureLocation[0,i] ) , \
                                        (fractureLocation[1,i] - fractureEdgeLength[0,i]/2) ,\
                                        (fractureLocation[2,i] + fractureEdgeLength[1,i]/2) ])
                
        # normal vector y-direction
        elif fractureOrientation[i] == 2:
            edgeLengthFactor = array([1 , 0 , 1]) # Does not seem to be used anywhere.
            fractureCornerPoints[i,0,:] = array([ (fractureLocation[0,i] + fractureEdgeLength[0,i]/2 ) , \
                                            (fractureLocation[1,i]) , \
                                            (fractureLocation[2,i] + fractureEdgeLength[1,i]/2 ) ])
            fractureCornerPoints[i,1,:] = array([ (fractureLocation[0,i] + fractureEdgeLength[0,i]/2 ) , \
                                            (fractureLocation[1,i]) , \
                                            (fractureLocation[2,i] - fractureEdgeLength[1,i]/2 ) ])
            fractureCornerPoints[i,2,:] = array([ (fractureLocation[0,i] - fractureEdgeLength[0,i]/2 ) , \
                                            (fractureLocation[1,i]) , \
                                            (fractureLocation[2,i] - fractureEdgeLength[1,i]/2 ) ])
            fractureCornerPoints[i,3,:] = array([ (fractureLocation[0,i] - fractureEdgeLength[0,i]/2 ) , \
                                            (fractureLocation[1,i]) , \
                                            (fractureLocation[2,i] + fractureEdgeLength[1,i]/2 ) ])
            
        # normal vector z-direction
        elif fractureOrientation[i] == 3:
            edgeLengthFactor = array([1 , 1 , 0]) # Does not seem to be used anywhere.
            fractureCornerPoints[i,0,:] = array([ (fractureLocation[0,i] + fractureEdgeLength[0,i]/2 ) , \
                                                  (fractureLocation[1,i] + fractureEdgeLength[1,i]/2 ) , \
                                                  (fractureLocation[2,i]) ]) 
            fractureCornerPoints[i,1,:] = array([ (fractureLocation[0,i] + fractureEdgeLength[0,i]/2 ) , \
                                                  (fractureLocation[1,i] - fractureEdgeLength[1,i]/2 ) , \
                                                  (fractureLocation[2,i]) ])
            fractureCornerPoints[i,2,:] = array([ (fractureLocation[0,i] - fractureEdgeLength[0,i]/2 ) , \
                                                  (fractureLocation[1,i] - fractureEdgeLength[1,i]/2 ) , \
                                                  (fractureLocation[2,i]) ])
            fractureCornerPoints[i,3,:] = array([ (fractureLocation[0,i] - fractureEdgeLength[0,i]/2 ) , \
                                                  (fractureLocation[1,i] + fractureEdgeLength[1,i]/2 ) , \
                                                  (fractureLocation[2,i]) ])
        
    # fracture size sorting for coloring of fractures
    fractureSize = fractureEdgeLength[0,:]*fractureEdgeLength[1,:]
       
    # round fracturePoints for discrete application of fractures
    fractureCornerPoints = around(fractureCornerPoints)
    # Here the coordinates of the fracture's corners are rounded to the nearest integer. 

    # cut the fractures at the boundaries
    fractureCornerPoints[fractureCornerPoints>xSize] = xSize
    fractureCornerPoints[fractureCornerPoints< 0.0] =  0.0
    
    if savePNG == 1: 
        # sort fracture sizes
        newRow = linspace(1,numberOfFractures,numberOfFractures)
        fractureSize = vstack([fractureSize, newRow])
        val = argsort(fractureSize[0,:]) 
        ind=fractureSize[:,val]
        ind = vstack([ind, newRow])
        val2 = argsort(ind[1,:]) 
        fractureSize=ind[:,val2]
           
        # color scheme for all fractures 
        mm = cm.ScalarMappable()
        ccc = mm.to_rgba(fractureSize[0,:]) # color fractures according to their sizes
            
        # plot fracture network
        dfn =  plt.figure()
        mng = plt.get_current_fig_manager()
        mng.resize(*mng.window.maxsize())    # have the figure use the entire screen.
        plt.hold(True)
        ax = Axes3D(dfn)
        ax.hold(True)
        plotPoints = 10


                
        # plot fracture network
        v = []
        for i in range(0,numberOfFractures):
        # Each fracture is plotted here using the four corner points and given a color corresponding to its size.
            x = [fractureCornerPoints[i,0,0], fractureCornerPoints[i,1,0], fractureCornerPoints[i,2,0], fractureCornerPoints[i,3,0]]
            y = [fractureCornerPoints[i,0,1], fractureCornerPoints[i,1,1], fractureCornerPoints[i,2,1], fractureCornerPoints[i,3,1]]
            z = [fractureCornerPoints[i,0,2], fractureCornerPoints[i,1,2], fractureCornerPoints[i,2,2], fractureCornerPoints[i,3,2]]
            v.append(zip(x, y, z))
            
        poly3dCollection = Poly3DCollection(v)
        poly3dCollection.set_color(ccc)
        poly3dCollection.set_edgecolor('k')
        sc = ax.add_collection3d(poly3dCollection)   
                                           
        # Here the plot layout, titles and such are set to ensure all resulting figures look the same.
        ax.set_xlabel('dimensionless x [-]')
        ax.set_ylabel('dimensionless y [-]')
        ax.set_zlabel('dimensionless z [-]')
        title_str = 'Fracture Network with %d fractures' % (numberOfFractures)
        ax.set_title(title_str)    
        x = 1; y = 1; width = 1200; height = 800    
        plt.grid()
        m = cm.ScalarMappable(cmap=cm.jet)
        m.set_array(fractureSize[0,:])
        cb = dfn.colorbar(m)
        plt.axis('equal')
        plt.axis('tight')
        ax.view_init(elev=18., azim=-43.)
        ax.set_xlim([xcoords[0]-round(fractureSizeMax/2), xcoords[1]+round(fractureSizeMax/2)])
        ax.set_ylim([ycoords[0]-round(fractureSizeMax/2), ycoords[1]+round(fractureSizeMax/2)])
        ax.set_zlim([zcoords[0]-round(fractureSizeMax/2), zcoords[1]+round(fractureSizeMax/2)])
        ax.hold(False)
        plt.hold(False)


        ## SAVE DATA AND FIGURES for GEOS and Graph Theory applications                
        saveFig = name+ '.png' # path and name of file we will save to.
        plt.savefig(workingDir +  saveFig)


    
    ## Adjust boundaries to include the inlet and outlet fractures in the GEOS simulation
    nx = 54    
    xcoords = array([-0.40, xSize+0.40])
        

    ############################################################################################
    ### removing fractures that are completely contained within other fractures
    ### and adding the full boundary.
    ############################################################################################

    redo = 1
    addBoundary = 1
    workingDir = workingDir + "plot/" 
    inputFile = name +'.xml'
    #def overlapCheck(workingDir, inputFile, redo, simulation, numberOfFracturesNeeded, addBoundary):
    # xmldoc     = minidom.parse( workingDir + inputFile ) 
    # itemlist = xmldoc.getElementsByTagName('Nodeset')
    # Inside of the xml file the nodesets are given as follows:
    # <Nodeset name="inlet" type="0"     
    #          xmin="-0.01 -0.01 1.99"     
    #          xmax="0.01 2.01 2.01" >
    # With these 3 pairs of coordinates all 8 cornerpoints of a fracture are described. 
    # The code above thus extracts these fracture coordinates from the xml file.
    # These commands initialize lists and the nBoxes counter which will be filled with the relevant information
    nBoxes= 0
    boxes    = []
    boxNames = []
   
    for i in range(0,numberOfFractures):
            # Here we run through each fracture (nodeset) in itemlist.
            x = fractureCornerPoints[i,2,0] - diff          # minimum coordinates
            y = fractureCornerPoints[i,2,1] - diff 
            z = fractureCornerPoints[i,2,2] - diff 
            minP  = box3D.PointT(x,y,z)                     # box3d.PointT turns these into an objects containing the 3 floats.
            x = fractureCornerPoints[i,0,0] + diff          # same for maximum coordinates
            y = fractureCornerPoints[i,0,1] + diff 
            z = fractureCornerPoints[i,0,2] + diff
            maxP  = box3D.PointT(x,y,z)
            box   = box3D.BoxT(minP, maxP, nBoxes)          # Using the minimum and maximum coordinates (minP and maxP) each fracture is represented by a box using box3D.BoxT.
            boxes.append( box )                             # This box is added to the boxes list
            nBoxes  += 1


    
    

    ############################################################################################
    # checking remaining fractures for those completely inside the space used by multiple other fractures
    ############################################################################################
    

    ############################################################################################
    ### Segmenting the entire domain using the outer coordinates of the fractures.
    ### Many of these segments will not correspond to a fracture.
    ############################################################################################
    intCoordx = []  # holds the minimum and maximum x-coordinates of the intersects
    intCoordy = []  # holds the minimum and maximum y-coordinates of the intersects
    intCoordz = []  # holds the minimum and maximum z-coordinates of the intersects
    fracID    = []  # holds the box number for the current coordinate

    for box in boxes:
       intCoordx.append(box.minP.x)
       intCoordx.append(box.maxP.x)
       intCoordy.append(box.minP.y)
       intCoordy.append(box.maxP.y)
       intCoordz.append(box.minP.z)
       intCoordz.append(box.maxP.z)
    
    # Reduce intersect coordinate lists to unique coordinate values    
    xCoord = np.unique(intCoordx)  # 
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
    for i in range(nBoxes):                                # Go through all boxes

        # X-coordinate intervals
        indLow = np.where(xCoord == intCoordx[2*i])         # Find the index of the first interval using xMin
        indUp  = np.where(xCoord == intCoordx[2*i+1])       # Find the index of the last  interval using xMax
        for j in range(int(indLow[0]), int(indUp[0])):     # Add the current fracture number to all intervals between xMin and xMax
            xID[j].append(i) 

        # Y-coordinate intervals
        indLow = np.where(yCoord == intCoordy[2*i])         # Find the index of the first interval using yMin
        indUp  = np.where(yCoord == intCoordy[2*i+1])       # Find the index of the last  interval using yMax
        for j in range(int(indLow[0]), int(indUp[0])):     # Add the current fracture number to all intervals between yMin and yMax
            yID[j].append(i)

        # Z-coordinate intervals
        indLow = np.where(zCoord == intCoordz[2*i])         # Find the index of the first interval using yMin
        indUp  = np.where(zCoord == intCoordz[2*i+1])       # Find the index of the last  interval using yMax
        for j in range(int(indLow[0]), int(indUp[0])):     # Add the current fracture number to all intervals between yMin and yMax
            zID[j].append(i)

    # find the segments that lie on fractures and turn them into nodes.
    nDomainSegBoxes = 0
    domainSegBoxes = []
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
                        minP  = box3D.PointT(xCoord[i]   , yCoord[j]   , zCoord[k])         # box3d.PointT turns these into an objects containing the 3 floats.
                        maxP  = box3D.PointT(xCoord[i+1] , yCoord[j+1] , zCoord[k+1])       # box3d.PointT turns these into an objects containing the 3 floats.
                        box   = box3D.BoxT(minP, maxP, nBoxes)                              # Using the minimum and maximum coordinates (minP and maxP) each segment is represented by a box using box3D.BoxT.
                        domainSegBoxes.append( box )                                        # This box is added to the boxes list
                        segments[k][j][i] =nDomainSegBoxes
                        #vI = g.add_vertex()
                        #cent[ vI ]  = [box.centroid().x, box.centroid().y, box.centroid().z]
                        if i == 0:
                            vertsInBox[0].append(nDomainSegBoxes)
                        elif i == len(xID)-1:
                            vertsInBox[1].append(nDomainSegBoxes)
                        nDomainSegBoxes +=1    


    ############################################################################################
    ### Checking which domain segments lie within fractures 
    ############################################################################################
    segmentsInBox = [] # lists which segments lie on a fracture
    segmentsInBoxCheck = [] # lists which segments lie on a fracture
    boxesInSeg    = [] # lists which fractures a segment belongs to
    boxesInSegCheck = []
    for b in boxes:
        segmentsInBox.append( [] )    #init
        segmentsInBoxCheck.append( [] )    #init
    for s in domainSegBoxes:
        boxesInSeg.append( [] )
        boxesInSegCheck.append( [] )

    for i in range(0, nDomainSegBoxes):
        for j in range(0, nBoxes):
                within = domainSegBoxes[i].inside(boxes[j]) # calculates coordinates of theoretic intersections.        
                if within == 1: # If this segments lies inside this fracture box
                    segmentsInBox[j].append( i )
                    boxesInSeg[i].append( j )

    # The logic of this part goes as follows:
    # If a fracture lies completely within space occupied by other fractures,
    # this fracture will not have a single segment to itsself.
    # All fractures that have at least 1 segment to themselves will be kept.
    # All others are removed.
    foundMultiple = []

    for i in range(0, nDomainSegBoxes):        
        if len(boxesInSeg[i])==1:
            foundMultiple.append(boxesInSeg[i])
            
    foundMultiple = np.unique(np.array(foundMultiple))

    remove = []
    for i in range(0,nBoxes):
        if i not in foundMultiple:
            remove.append(i)


    found = len(remove)
    

    



    ############################################################################################
    # Adding domain boundary fractures
    ############################################################################################
    if addBoundary == 1:

        leftBC  = 0
        rightBC = 0
        nfound  = 0
        BC = []
        used = 0
        for i, box in enumerate(boxes):
            if i not in remove and used < numberOfFracturesNeeded:
                used += 1
                if used == numberOfFracturesNeeded:
                    break
                elif box.minP.x < 0.0 and box.maxP.x > 0.02:
                    minP = box3D.PointT(-0.41, box.minP.y, box.minP.z)
                    maxP = box3D.PointT( 0.01, box.maxP.y, box.maxP.z)
                    newBox = box3D.BoxT(minP, maxP, nfound)
                    BC.append(newBox)
                    leftBC += 1
                    nfound += 1
                elif (box.minP.x < (xSize -2*diff)) and (box.maxP.x > xSize):
                    minP = box3D.PointT( xSize-diff, box.minP.y, box.minP.z)
                    maxP = box3D.PointT(xSize+0.4+diff, box.maxP.y, box.maxP.z)
                    newBox = box3D.BoxT(minP, maxP, nfound)
                    BC.append(newBox)
                    rightBC += 1
                    nfound  += 1
                

        print( "left ",used )
        assert numberOfFracturesNeeded == used


        
     


        ############################################################################################
        ### Segmenting of the entire domain using the coordinates of the boundary fractures 
        ### Many of these segments will not correspond to a fracture.
        ############################################################################################
        intCoordx = []  # holds the minimum and maximum x-coordinates of the intersects
        intCoordy = []  # holds the minimum and maximum y-coordinates of the intersects
        intCoordz = []  # holds the minimum and maximum z-coordinates of the intersects
        fracID    = []  # holds the box number for the current coordinate

        for box in BC:
           intCoordx.append(box.minP.x)
           intCoordx.append(box.maxP.x)
           intCoordy.append(box.minP.y)
           intCoordy.append(box.maxP.y)
           intCoordz.append(box.minP.z)
           intCoordz.append(box.maxP.z)
        
        # Reduce intersect coordinate lists to unique coordinate values    
        xCoord = np.unique(intCoordx)  # 
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
        for i in range(nfound):                                # Go through all boxes

            # X-coordinate intervals
            indLow = np.where(xCoord == intCoordx[2*i])         # Find the index of the first interval using xMin
            indUp  = np.where(xCoord == intCoordx[2*i+1])       # Find the index of the last  interval using xMax
            for j in range(int(indLow[0]), int(indUp[0])):     # Add the current fracture number to all intervals between xMin and xMax
                xID[j].append(i) 

            # Y-coordinate intervals
            indLow = np.where(yCoord == intCoordy[2*i])         # Find the index of the first interval using yMin
            indUp  = np.where(yCoord == intCoordy[2*i+1])       # Find the index of the last  interval using yMax
            for j in range(int(indLow[0]), int(indUp[0])):     # Add the current fracture number to all intervals between yMin and yMax
                yID[j].append(i)

            # Z-coordinate intervals
            indLow = np.where(zCoord == intCoordz[2*i])         # Find the index of the first interval using yMin
            indUp  = np.where(zCoord == intCoordz[2*i+1])       # Find the index of the last  interval using yMax
            for j in range(int(indLow[0]), int(indUp[0])):     # Add the current fracture number to all intervals between yMin and yMax
                zID[j].append(i)

        # find the segments that lie on fractures and turn them into nodes.
        nDomainSegBoxes = 0
        domainSegBoxes = []
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
                            minP  = box3D.PointT(xCoord[i]   , yCoord[j]   , zCoord[k])         # box3d.PointT turns these into an objects containing the 3 floats.
                            maxP  = box3D.PointT(xCoord[i+1] , yCoord[j+1] , zCoord[k+1])       # box3d.PointT turns these into an objects containing the 3 floats.
                            box   = box3D.BoxT(minP, maxP, nBoxes)                              # Using the minimum and maximum coordinates (minP and maxP) each segment is represented by a box using box3D.BoxT.
                            domainSegBoxes.append( box )                                        # This box is added to the boxes list
                            segments[k][j][i] =nDomainSegBoxes
                            #vI = g.add_vertex()
                            #cent[ vI ]  = [box.centroid().x, box.centroid().y, box.centroid().z]
                            if i == 0:
                                vertsInBox[0].append(nDomainSegBoxes)
                            elif i == len(xID)-1:
                                vertsInBox[1].append(nDomainSegBoxes)
                            nDomainSegBoxes +=1    


        ############################################################################################
        ### Checking which domain segments lie within fractures 
        ### and which domain segments additionally lie within an intersect.
        ############################################################################################
        segmentsInBox = [] # lists which segments lie on a fracture
        boxesInSeg    = [] # lists which fractures a segment belongs to
        for b in BC:
            segmentsInBox.append( [] )    #init
        for s in domainSegBoxes:
            boxesInSeg.append( [] )
    
        for i in range(0, nDomainSegBoxes):
            for j in range(0, nfound):
                    within = domainSegBoxes[i].inside(BC[j]) # calculates coordinates of theoretic intersections.        
                    if within == 1: # If this segments lies inside this fracture box
                        segmentsInBox[j].append( i )
                        boxesInSeg[i].append( j )

        foundMultiple = []
    

        for i in range(0, nDomainSegBoxes):
            
            if len(boxesInSeg[i])==1:
                foundMultiple.append(boxesInSeg[i])
                
        foundMultiple = np.unique(np.array(foundMultiple))


        removeBC = []
        for i in range(0,nfound):
            if i not in foundMultiple:
                removeBC.append(i)

        found = found + len(removeBC)


        if redo == 0:
            assert len(removeBC) == 0



    if redo == 1:      


        # read template
        if GEOSjob == 1:
            template = 'template_for_removal.xml'
        else:
            template = 'templateForGraph.xml'

        f        = open(template,'r'); # read in template
        lines    = f.readlines()
        f.close()

        filename = inputFile.split(".")[0] + '.xml'  
        print( workingDir + filename
        
        fileID = open(workingDir + filename,'w'); # open file for writing; append below existing contents
        for line in lines: # write template lines to new file
            fileID.write(line)
            
        # NODESETS
        # The final part saves the nodes of the fracture-plane segments. 
        # The 3 sets of coordinates can be used to make 8 cornerpoints of a thin rectangular box (a fracture). 

        # boundary fractures
        fileID.write('<!--boundary nodesets --> \n')
        for i, box in enumerate(BC):
            if i not in removeBC:
                if GEOSjob == 1:
                    fileID.write('    <Nodeset name="fracturePlane" type="0"     \n')
                    fileID.write('         xmin="%1.2f %1.2f %1.2f"        \n'% (box.minP.x, box.minP.y, box.minP.z) )
                else:
                    fileID.write('<Nodeset xmin="%1.2f %1.2f %1.2f"        \n'% (box.minP.x, box.minP.y, box.minP.z) )

                fileID.write('         xmax="%1.2f %1.2f %1.2f"/>      \n'% (box.maxP.x, box.maxP.y, box.maxP.z) )

        # fractures
        used = 0
        fileID.write('\n<!--normal nodesets --> \n')
        for i, box in enumerate(boxes):
            if i not in remove and used < numberOfFracturesNeeded:
                if GEOSjob == 1:
                    fileID.write('    <Nodeset name="fracturePlane" type="0"     \n')
                    fileID.write('         xmin="%1.2f %1.2f %1.2f"        \n'% (box.minP.x, box.minP.y, box.minP.z) )
                else:
                    fileID.write('<Nodeset xmin="%1.2f %1.2f %1.2f"        \n'% (box.minP.x, box.minP.y, box.minP.z) )

                fileID.write('         xmax="%1.2f %1.2f %1.2f"/>      \n'% (box.maxP.x, box.maxP.y, box.maxP.z) )
                used += 1
   




        

        fileID.write('</Nodesets>\n');

        if GEOSjob == 1:
            Pin  = 2.0e6           # inlet pressure
            Pout = 1.0e6           # outlet pressure
            fileID.write('\n')
            fileID.write('<BoundaryConditions> \n')
            fileID.write('    <BoundaryCondition  object="Node" fieldname="Displacement" setnames="xneg" component="0" scale="0.0" fieldtype="Vector" />  \n')
            fileID.write('    <BoundaryCondition  object="Node" fieldname="Displacement" setnames="xpos" component="0" scale="0.0" fieldtype="Vector" /> \n')
            fileID.write('    <BoundaryCondition  object="Node" fieldname="Displacement" setnames="yneg" component="1" scale="0.0" fieldtype="Vector" /> \n')
            fileID.write('    <BoundaryCondition  object="Node" fieldname="Displacement" setnames="ypos" component="1" scale="0.0" fieldtype="Vector" /> \n')
            fileID.write('    <BoundaryCondition  object="Node" fieldname="Displacement" setnames="zneg" component="2" scale="0.0" fieldtype="Vector" />  \n')
            fileID.write('    <BoundaryCondition  object="Node" fieldname="Displacement" setnames="zpos" component="2" scale="0.0" fieldtype="Vector" /> \n')
            fileID.write('<!-- \n')
            fileID.write('FLUID FLOW BC \n')
            fileID.write('--> \n')
            fileID.write('    <!--<BoundaryCondition object="Face" fieldname="VolumeRate" setnames="inlet" scale="1.0e-9" table="ttable" option="1" /> /--> \n')
            fileID.write('    <BoundaryCondition object="Face" fieldname="Pressure" setnames="inlet"  scale="{0}" timetable="ttable" />   \n'.format(Pin))
            fileID.write('    <BoundaryCondition object="Face" fieldname="Pressure" setnames="outlet" scale="{0}" timetable="ttable" />   \n'.format(Pout))
            fileID.write('</BoundaryConditions> \n')
            fileID.write('\n')
 
 
            fileID.write('<Output writePlot="1"  \n')
            fileID.write('  plot_interval="$:plotInterval" \n')
            fileID.write('  writeRestart="0"  \n')
            fileID.write('  restart_interval="1000.0"  \n')
            fileID.write('  plotfile_root="plot{0}"  \n'.format(simulation))
            fileID.write('  parallel_silo="16"  \n')
            fileID.write('  slave_directory="sub" \n')
            fileID.write('  writeFEMFaces="1"  \n')
            fileID.write('  writeFEMEdges="0" /> \n')
            fileID.write('\n')
            fileID.write('</Problem>')

        fileID.close()
