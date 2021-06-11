#! /usr/bin/python

# Copyright:
# 	Alex Hobe
# 	Daniel Vogler
# 	Martin P. Seybold


# TODO
# - extract abstractions?
#	- import abstractions
# - apeture value needed throughout.
# 	- constant now, readable from file in a non-existent future.
# - produce inheritance?

import sys # import sys.exit only?
from xml.dom import minidom
import numpy as np

def abstractFracGeomGetter( fileName,fileType ):
	""" Set of input methods, that extract fracture geometries from a file and convert this to the pipeline assumptions/abstraction.

	"""

def getFractureBoxes( fileName, fileType ):
	# fileType as input, as multiple csv and xml types possible.
	# use inheritence instead of hardcoded if/elif
	if fileType == "GEOSxml":
		#fracBoxes = readGEOSxml( fileName )
		intCoordx, intCoordy, intCoordz, nBoxes = readGEOSxml( fileName )
	elif fileType == "fracCSV":
		print( "fracCSV reader not yet implemented" )
		# TODO: error/exception instead of print
		sys.exit()
		fracBoxes = readFracCsv( fileName )
	
	# return fracBoxes
	return intCoordx, intCoordy, intCoordz, nBoxes

def abstractFracGeomReader( fileName ):
	""" Set of readers, that extract fracture geometries from a file.

		Usage: Input fileID, output fracture geometries
	"""

def readGEOSxml( fileName ):
	coordList = getCoordListGEOSxml( fileName )
	#fracBoxes = makeFracBoxGEOSxml( coordList )
	intCoordx, intCoordy, intCoordz, nBoxes = makeFracBoxGEOSxml( coordList )
	# return fracBoxes
	return intCoordx, intCoordy, intCoordz, nBoxes
	
def abstractCoordListGetter( fileName ):
	""" Set of converters between fracture-geometry-file specifics and fractureBox implementation.
	"""

def getCoordListGEOSxml( fileName ):
		xmldoc = minidom.parse( fileName ) 
		coordList = xmldoc.getElementsByTagName( 'Nodeset' )
		return coordList


def abstractFracBoxMaker( coordList ):
	""" Set of converters between fracture-geometry-file specifics and fractureBox implementation.
	"""

def makeFracBoxGEOSxml( coordList ):
	# TODO: 
	# - adjust to using box3D < not optimized for the problem?
	# - provide aperture value elsewhere 
	# 	- first needed here. Export as value of fracBoxes object?
	# - adjust to fracBoxes first
	# 	- fracBoxes.intCoordx etc.

	# initialize coord minMax lists
	intCoordx = []    
	intCoordy = []    
	intCoordz = []
	aperture   = 1.0e-5
	diff       = 0.01 - aperture # used to adjust for the arbitrary overlap to capture the nodes in GEOS.
	nBoxes   = 0

	for s in coordList:
		x,y,z = s.attributes['xmin'].value.split()      # The minimum coordinates are extracted and
		intCoordx.append(float(x) + diff)				# added to the coordinate lists
		intCoordy.append(float(y) + diff)
		intCoordz.append(float(z) + diff)
		x,y,z = s.attributes['xmax'].value.split()      # same for maximum coordinates
		intCoordx.append(float(x) - diff)
		intCoordy.append(float(y) - diff)
		intCoordz.append(float(z) - diff)
		nBoxes  += 1

	#return fracBoxes
	return intCoordx, intCoordy, intCoordz, nBoxes

def makeFracBoxGEOSxmlForISPM( coordList ):
	# This differs slightly from HSPM version.
	# Can we combine the two?
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
    
    for s in coordList:
    # Here we run through each fracture (nodeset) in coordList.
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
    return xCoord, yCoord, zCoord, nBoxes, longDimArray, verticesInBox, widthArray

# TODO:
# - extend getFractureBoxes using inheritance here.

def readFracCsv( abstractFracGeomReader ):
	# TODO:
	# - create desired csv
	# - implement reader
	coordList = getCoordListFracCsv( fileName )
	fracBoxes = makeFracBoxesFracCsv( coordList )
	return fracBoxes





