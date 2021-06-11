import unittest
import numpy as np
import sys
sys.path.insert(1, '../')
from graphMaking import *

if __name__ == '__main__':
	unittest.main()

class TestGraohMaking(unittest.TestCase):
	# given
	def setUp(self):
		self.intCoordx = [-0.40001, -0.19999, 10.19999, 10.40001]
		self.intCoordy = [-9.999999999999593e-06, 10.00001, -9.999999999999593e-06, 10.00001]
		self.intCoordz = [-9.999999999999593e-06, 10.00001, -9.999999999999593e-06, 10.00001]
		self.nBoxes = 2
	def test_segmentDomain(self):
		# when
		xID, yID, zID, xCoord, yCoord, zCoord = segmentDomain( self.intCoordx, self.intCoordy, self.intCoordz, self.nBoxes )
		
		# then
		xIDcheck = [[], [], []]
		yIDcheck = [[]]
		zIDcheck = [[]]
		xCheck = np.array([-0.40001, -0.19999, 10.19999, 10.40001])
		yCheck = np.array([-1.000000e-05,  1.000001e+01])
		zCheck = np.array([-1.000000e-05,  1.000001e+01])
		self.assertEqual( xID, xIDcheck )
		self.assertEqual( yID, yIDcheck )
		self.assertEqual( zID, zIDcheck )
		self.assertEqual( xCoord.all(), xCheck.all() )
		self.assertEqual( yCoord.all(), yCheck.all() )
		self.assertEqual( zCoord.all(), zCheck.all() )
