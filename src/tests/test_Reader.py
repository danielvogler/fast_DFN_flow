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
	
	def test_getFractureBoxesISPM(self):
		# when
		itemlist = getCoordListGEOSxml(self.geosXMLinputFile)
		xCoord, yCoord, zCoord, nBoxes, longDimArray, verticesInBox, widthArray = makeFracBoxGEOSxmlForISPM( itemlist )       
		
		# then
		xCheck = [[-0.41, 10.19], [-0.19, 10.41]]
		yCheck = [[-0.01, -0.01], [10.01, 10.01]]
		zCheck = [[-0.01, -0.01], [10.01, 10.01]]
		nCheck = 2
		longDimCheck = [1, 1] 
		verticesCheck = [[], []] 
		widthCheck = [10.02, 10.02]
		self.assertEqual( xCoord, xCheck )
		self.assertEqual( yCoord, yCheck )
		self.assertEqual( zCoord, zCheck )
		self.assertEqual( nBoxes, nCheck )
		self.assertEqual( longDimArray, longDimCheck )
		self.assertEqual( verticesInBox, verticesCheck )
		self.assertEqual( widthArray, widthCheck )

	# def test_getCoordListGEOSxml(self):
	#	# No idea how to test this, other than the previous test.
	#	# when
	#	coordList = getCoordListGEOSxml( self.geosXMLinputFile )

		





