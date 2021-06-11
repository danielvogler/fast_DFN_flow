import unittest
import sys
sys.path.insert(1, '../')
from getFracGeometries import *

if __name__ == '__main__':
	unittest.main()

class TestReaders(unittest.TestCase):
	# given
	def setUp(self):
		self.geosXMLinputFile = "testFile.xml"
		self.geosXMLType = "GEOSxml"
	def test_getFractureBoxes(self):
		# when
		intCoordx, intCoordy, intCoordz, nBoxes = getFractureBoxes(self.geosXMLinputFile, self.geosXMLType)
		
		# then
		xCheck = [-0.40001, -0.19999, 10.19999, 10.40001]
		yCheck = [-9.999999999999593e-06, 10.00001, -9.999999999999593e-06, 10.00001]
		zCheck = [-9.999999999999593e-06, 10.00001, -9.999999999999593e-06, 10.00001]
		nCheck = 2
		self.assertEqual( intCoordx, xCheck )
		self.assertEqual( intCoordy, yCheck )
		self.assertEqual( intCoordz, zCheck )
		self.assertEqual( nBoxes, nCheck )

	# def test_getCoordListGEOSxml(self):
	#	# No idea how to test this, other than the previous test.
	#	# when
	#	coordList = getCoordListGEOSxml( self.geosXMLinputFile )

		





