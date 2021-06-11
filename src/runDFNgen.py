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
import scipy.stats as stats

import box3D
import graph_tool.all as gt # https://graph-tool.skewed.de/static/doc/index.html
from overlap_and_boundary import *
from generateMsub import *
from DFNgenFunction import *
import os
import time


############################
### User input
############################

# Select range and number of simulation for the fracture density distribution (normal distribution between lower and upper cutoff values)
lower =  50
upper = 450
mu    = 250
sigma =  50
nIter = 1000
trymax = 100

# Restart generation from point in file.
Restart   = 0   # 0 for start a new, 1 for Restart from save point.
savePoint = 1000  # choose after how many iterations a savepoint occurs

# Create GEOS job files?
GEOSjob = 0  # 1 for create jobfiles, 0 for not creating them.
# Only gets done, if Restart == 0.

# Choose to time the generation process
genTime = 1 # 1 time generation, 0 leave out.

# Choose to save the DFN picture
savePNG = 0

############################
############################
############################




## Check which benchmark number we are at:
filename = "currentBenchmark.txt"
f        = open(filename,'r'); # read in benchmark number
nBenchmark = int(f.readline())
f.close()
print( "Creating files for benchmark: ", nBenchmark )





## Create folders and filenames from input parameters.    
large = nIter - 1
numberOfFracturesNeeded = np.round( (upper + lower)/2 )
print( numberOfFracturesNeeded )
simName = "bm" + str(nBenchmark) + "_f" + str(numberOfFracturesNeeded)
path_fig = "../../figures"
workingDir = "../benchmark/" + simName + "/"

if not os.path.exists(workingDir):
          os.makedirs(workingDir)

if not os.path.exists(workingDir + "plot"):
          os.makedirs(workingDir + "plot")

## Check if new generation or restarting one.
if Restart == 1:
    filename = workingDir + "Restart.txt"
    f        = open(filename,'r'); # read in benchmark number
    small = int(f.readline())
    f.close()
    print( "Restarting from simulation nr: ", small  )
else:
    small = 0
    
    if GEOSjob == 1:
        # create bash script to submit all runs
        SUBMIT(workingDir, simName, small, large)

##########
# used to only generate the bash scripts in case something needs to be changed there...
##########
#for simulation in range(small, large+1):
#    # create bash scripts to run the simulations in GEOS
#    if GEOSjob == 1:
#        MSUB(workingDir, simName + "_{0}".format(simulation), simulation)
#        SUBMIT(workingDir, simName, small, large)
# a = 1/0
##########

# Read in previously created fracture density distribution or create a new one using the input parameters.
if Restart == 1: 
    filename = workingDir + "rhoFrac.txt"    
    fileID = open(filename,'r'); # open file for reading
    rhoFracArray =  np.genfromtxt (fileID, delimiter=",")
else:
    X = stats.truncnorm(
    (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
    rhoFracArray = X.rvs(nIter)
    assert np.min(rhoFracArray)>0
    
    for rho in rhoFracArray:
        rho = int(np.round(rho))


    filename = workingDir + "rhoFrac.txt"    
    fileID = open(filename,'w'); # open file for writing; discard existing contents
    for i in range(small,large+1):
        if i == 0:
            fileID.write("{0}".format(rhoFracArray[i]))
        else:
            fileID.write(",{0}".format(rhoFracArray[i]))


# Initialize file to time the generation process.
if genTime == 1:
    filename = workingDir + "generation.txt"
    f        = open(filename,'w'); # read in benchmark number
    f.write('nFrac,t[s]\n')
    f.close()

# Initialize file that takes care of failed generations
filename = workingDir + "failedDFN.txt"
f        = open(filename,'w') # read in benchmark number
f.write('nDFN,nFrac\n')
f.close()


##############################################################
### Loop through fracture density distribution and create DFN.
##############################################################

for simulation in range(small, large+1):
    t = time.time()
    numberOfFracturesNeeded = int(rhoFracArray[simulation])
    #numberOfFracturesNeeded = 367

    # Generate the DFN with the desired fracture density.
    done = 0
    for tryout in range(0,trymax):
        print( "trying ", tryout )
        if done == 0:
            try:
                generateDFN(simulation, workingDir, simName, numberOfFracturesNeeded, savePNG, GEOSjob)
                done = 1
                break

            except Exception as e:
                if tryout == trymax-1:
                    print( "Failed {0} times.".format(trymax) )
                    # Save the iteration that failed to generate
                    filename = workingDir + "failedDFN.txt"
                    f        = open(filename,'a') # read in benchmark number
                    f.write('{0}'.format(simulation))
                    f.write(',{0}\n'.format(numberOfFracturesNeeded))
                    f.close()


    elapsed = time.time() - t
    print( numberOfFracturesNeeded, elapsed )
    print( "" )
    

    # create bash scripts to run the simulations in GEOS
    if GEOSjob == 1:
        MSUB(workingDir, simName + "_{0}".format(simulation), simulation)

    # save elapsed time to file
    if genTime == 1:
        filename = workingDir + "generation.txt"
        f        = open(filename,'a') # read in benchmark number
        f.write('{0}'.format(numberOfFracturesNeeded))
        f.write(',{0}\n'.format(elapsed))
        f.close()

    # create savePoint to start from if the generation crashes completely.
    if simulation%savePoint==0:
        filename = workingDir + "Restart.txt"
        f        = open(filename,'w'); # read in benchmark number
        f.write('{0}'.format(simulation))
        f.close()

sys.exit() 
# If all DFN were created successfully, update the benchmark counter.
nBenchmark += 1
filename = "currentBenchmark.txt"
f        = open(filename,'w'); # overwrite file to input new benchmark number
f.write('{0}'.format(nBenchmark))
f.close()
print( "Next benchmark will be nr ", nBenchmark )

