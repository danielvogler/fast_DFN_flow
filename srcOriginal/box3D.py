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

__author__ = "martin"
__date__ = "$Jun 24, 2015 1:03:06 PM$"

if __name__ == "__main__":
    print "Just executed box3D.py"

class PointT(object): # Creates a point object with three float coordinates
    #@profile
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x=float(x)
        self.y=float(y)
        self.z=float(z)
    #@profile
    def str(self):
        return( "(" + str(self.x) +", "+ str(self.y) +", "+ str(self.z )+")")
     
class BoxT(object): # Creates a box object 
    #@profile
    def __init__(self, minP, maxP, idx): # The 6 coordinates in minP and maxP can be used to determine the 8 corner points
        self.minP = minP
        self.maxP = maxP
        self.idx  = idx  # gives the box a number
    #@profile
    def isEmpty(self):   # checks if there is a negative size of the box (Returns True). Returns False if all dimensions are positive.
        return( (self.minP.x > self.maxP.x) or (self.minP.y > self.maxP.y) or (self.minP.z > self.maxP.z) ) 
    #@profile
    def intersect(self, other): # defines the 6 coordinates of the intersect box between two fractures with which the 8 corner points can be found.      
        minP = PointT( max(self.minP.x, other.minP.x),
                       max(self.minP.y, other.minP.y),
                       max(self.minP.z, other.minP.z)  )
        maxP = PointT( min(self.maxP.x, other.maxP.x),
                       min(self.maxP.y, other.maxP.y),
                       min(self.maxP.z, other.maxP.z)  )           
        return( BoxT(minP, maxP, self.idx) ) # returns the intersect as a box.

    #@profile
    def intersectSegment(self, other): # defines the 6 coordinates of the intersect box between two fractures with which the 8 corner points can be found.
        eps      = 1e-13 # used to check for machine precision errors.
        if  abs(self.minP.x - other.minP.x)<= eps and \
            abs(self.maxP.x - other.maxP.x)<= eps and \
            abs(self.minP.y - other.minP.y)<= eps and \
            abs(self.maxP.y - other.maxP.y)<= eps and \
            abs(self.minP.z - other.minP.z)<= eps and \
            abs(self.maxP.z - other.maxP.z)<= eps:
            print "equal fracture found"
            minP = PointT(2,2,2)
            maxP = PointT(1,1,1)
            equal = 1               # This variable is used to later ignore all intersects with this fracture copy.
        else:
            minP = PointT( max(self.minP.x, other.minP.x),
                           max(self.minP.y, other.minP.y),
                           max(self.minP.z, other.minP.z)  )
            maxP = PointT( min(self.maxP.x, other.maxP.x),
                           min(self.maxP.y, other.maxP.y),
                           min(self.maxP.z, other.maxP.z)  )
            equal = 0

           
        return( BoxT(minP, maxP, self.idx), equal ) # returns the intersect as a box.

    def equalCheck(self, other):  # checks if two fractures are equal to each other
        eps      = 1e-13 # used to check for machine precision errors.
        if  abs(self.minP.x - other.minP.x)<= eps and \
            abs(self.maxP.x - other.maxP.x)<= eps and \
            abs(self.minP.y - other.minP.y)<= eps and \
            abs(self.maxP.y - other.maxP.y)<= eps and \
            abs(self.minP.z - other.minP.z)<= eps and \
            abs(self.maxP.z - other.maxP.z)<= eps:
            equal = 1               # This variable is used to later ignore all intersects with this fracture copy.
        else:
            equal = 0
        return equal

    #@profile
    def inside(self,other): # checks if one box is found completely within another box.
        eps      = 1e-5 # used to check for machine precision errors.        
        if  other.minP.x - eps <= self.minP.x and \
            other.maxP.x + eps >= self.maxP.x and \
            other.minP.y - eps <= self.minP.y and \
            other.maxP.y + eps >= self.maxP.y and \
            other.minP.z - eps <= self.minP.z and \
            other.maxP.z + eps >= self.maxP.z:
            ok = 1
            return ok
        else:
            ok = 0
            return ok

    #@profile
    def centroid(self): # calculates the geometric center
        return( PointT( 0.5*(float(self.minP.x) + float(self.maxP.x) ),
                        0.5*(float(self.minP.y) + float(self.maxP.y) ),
                        0.5*(float(self.minP.z) + float(self.maxP.z) ) ) )
    #@profile
    def projectLongestDim(self, point ): # searches for the longest dimension of a box and returns the corresponding coordinate value.
        dx , dy, dz = self.spread()
        if   dx >= dy and dx >= dz:
            return( point.x )
        elif dy >= dx and dy >= dz:
            return( point.y )
        else:
            return( point.z )

    #@profile
    def projectMediumDim(self, point ): # searches for the longest dimension of a box and returns the corresponding coordinate value.
        dx , dy, dz = self.spread()
        if   dx >= dy and dx <= dz:
            return( point.x )
        elif dy >= dx and dy <= dz:
            return( point.y )
        else:
            return( point.z )
    #@profile
    def ShortestDim(self ): # searches for the shortest dimension of a box and returns the corresponding coordinate value.
        dx , dy, dz = self.spread()
        if   dx <= dy and dx <= dz:
            return( 0 )     # will access a list with this index
        elif dy <= dx and dy <= dz:
            return( 1 )
        else:
            return( 2 )
    #@profile
    def LongestDim(self ): # searches for the shortest dimension of a box and returns the corresponding coordinate value.
        dx , dy, dz = self.spread()
        if   dx >= dy and dx >= dz:
            return( 0 )
        elif dy >= dx and dy >= dz:
            return( 1 )
        else:
            return( 2 )
    #@profile
    def MediumDim(self ): # searches for the shortest dimension of a box and returns the corresponding coordinate value.
        short = self.ShortestDim()
        Long  = self.LongestDim()
        if short != 0 and Long != 0:
            return( 0 )
        elif short != 1 and Long != 1:
            return( 1 )
        else:
            return( 2 )
                        
    #@profile
    def spread(self): # calculates the dimensions of the box in x, y, and z-direction.
        return( [self.maxP.x - self.minP.x,
                 self.maxP.y - self.minP.y,
                 self.maxP.z - self.minP.z ] )
    #@profile
    def length(self): # defined as the longest dimension
        return( sorted( self.spread() )[2] )
    #@profile
    def width(self):  # defined as the intermediate dimension
        return( sorted( self.spread() )[1] )
    #@profile
    def height(self): # defined as the shortest dimension (usually the aperture).
        return( sorted( self.spread() )[0])
    #@profile    
    def maxSpread(self): # returns the longest dimension of a box.
        return( max( float(self.maxP.x) - float(self.minP.x),
                     float(self.maxP.y) - float(self.minP.y),
                     float(self.maxP.z) - float(self.minP.z) ) )
    #@profile
    def str(self): # returns the important data associated with a box as a string.
        return("box.idx=" + str(self.idx)
            + "\tminP" + self.minP.str()
            + "\tmaxP" + self.maxP.str() )

        
        
        

