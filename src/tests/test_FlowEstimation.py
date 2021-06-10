import unittest
import sys
sys.path.insert(1, '../')
from HSPM_opt import segmentation
from ISPM_opt import simpleMethod

class Test2Dcase(unittest.TestCase):
	# better to use a separate file specifically for testing

	# given
	def setUp(self):
		self.workingDir = "../../benchmark/bm54_f17/plot/"
		self.inputFile  =              "bm54_f17_3.xml"
	def test_HSPM_opt_2D(self):
		# when
		Q, numPaths, nVertices, nEdges = segmentation(self.workingDir, self.inputFile)

		# then
		Q2D =  2.159524473137551e-08
		numPaths2D = 14
		nVertices2D = 259
		nEdges2D = 610	
		self.assertEqual(numPaths, numPaths2D)
		self.assertEqual(nVertices, nVertices2D)
		self.assertEqual(nEdges, nEdges2D)
		self.assertEqual(Q, Q2D) # not sure why this float comparison works.

	def test_ISPM_opt_2D(self):
		# when
		Q, numPaths, nVertices, nEdges = simpleMethod(self.workingDir, self.inputFile)

		# then
		Q2D =  1.0980362356816208e-08
		numPaths2D = 5
		nVertices2D = 26
		nEdges2D = 53	
		self.assertEqual(numPaths, numPaths2D)
		self.assertEqual(nVertices, nVertices2D)
		self.assertEqual(nEdges, nEdges2D)
		self.assertEqual(Q, Q2D) # not sure why this float comparison works.



class Test3Dcase(unittest.TestCase):
	# better to use a separate file specifically for testing

	# given
	def setUp(self):
		self.workingDir = "../../benchmark/bm95_f250/plot/"
		self.inputFile  =              "bm95_f250_0_2.xml"

	def test_HSPM_opt_3D(self):
		# when
		Q, numPaths, nVertices, nEdges = segmentation(self.workingDir, self.inputFile)

		# then
		Q3D =  3.0162195648101415e-07
		numPaths3D = 40
		nVertices3D = 4800
		nEdges3D = 10314	
		self.assertEqual(numPaths, numPaths3D)
		self.assertEqual(nVertices, nVertices3D)
		self.assertEqual(nEdges, nEdges3D)
		self.assertEqual(Q, Q3D) # not sure why this float comparison works.

	def test_ISPM_opt_3D(self):
		# when
		Q, numPaths, nVertices, nEdges = simpleMethod(self.workingDir, self.inputFile)

		# then
		Q3D =  3.060978751842549e-07
		numPaths3D = 178
		nVertices3D = 1495
		nEdges3D = 5386	
		self.assertEqual(numPaths, numPaths3D)
		self.assertEqual(nVertices, nVertices3D)
		self.assertEqual(nEdges, nEdges3D)
		self.assertEqual(Q, Q3D) # not sure why this float comparison works.



		
#		Q, numPaths, nVertices, nEdges = segmentation(workingDir, inputFile)

if __name__ == '__main__':
	unittest.main()
