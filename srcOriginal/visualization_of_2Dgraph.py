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
__date__ = "$Jun 21, 2015 2:57:02 PM$"

import numpy as np
import graph_tool.all as gt     # import networkx as nx  # has mincost flow
import graph_tool.topology as tp
import sys
import copy as cp
import os
import csv

normalize = 0  # print normalized result = 1, print actual result = 0.
get_path  = 1  # Extract path information

eps  = 1.0e-13 # used to check for machine precision errors.
diff = 0.02    # used to adjust the arbitrary overlap to capture the nodes in GEOS.


#######
### Choose which version of the algorithm
#######
#version = "per"
#version = "width"
version = "seg"
#version = "seg002"


#######
### Choose simulation folder and name
#######
name    = "bm54_f17"
simName = name + "_3"
workingDir      = "../benchmark/"   +    name + "/plot/"
inputFile       = workingDir        + simName +"_Graph_" + version + ".xml.gz"
outputPicGraph  = ".." + inputFile.split(".")[2] + "_checking.png"
outputpicFlow   = ".." + inputFile.split(".")[2] + "_Flow.png"
outputpicMCFlow = ".." + inputFile.split(".")[2] + "_MC_Flow.png"

V_PressureFile = ""

g = gt.Graph()
g.load( inputFile ) # loads the graph created in geos2graph

pos       = g.vertex_properties["pos"]      # extracts node positions
vLabel    = g.vertex_properties["vLabel"]   # extracts node names (box 1 nr. & box 2 nr.)
length    = g.edge_properties["length"]     # extracts length of each fracture segment
width     = g.edge_properties["width"]      # extracts width of each fracture segment
eFromBox  = g.edge_properties["eFromBox"]   # extracts the number of the original fracture to which each segment belongs.
apertures = g.graph_properties["apertures"] # extracts the aperture of each fracture segment
cent      = g.vertex_properties["cent"]     # extracts the coordinates of the vertices for the tortuosity calculations.

# setup edge capacities for maxFlow
cap    = g.new_edge_property("double")      # capacity for flow in each edge
eLabel = g.new_edge_property("string")      # will be used to label each edge in the final figures.

minVal =  10000.0 
maxVal = -10000.0                           # Initializes these variables so that the smallest and largest edge length can be determined.
mu     = 0.001
cubicLawConst = 12
cap_max = 0     
stuff = 0                            # Find the maximum capacity to normalize

a=1.0e-5

for e in g.edges():
    minVal = min( minVal, length[e])        # determines the smallest edge in the graph
    maxVal = max( maxVal, length[e])        # determines the lagest edge in the graph (both do not seem to be used anywhere.)
    if length[e] <= 1e-9:                   # set capacity for overlaping nodes.
        cap[e] = 100.0
    else :
        cap[e] = float(width[e] * a**3 / (length[e] * mu* cubicLawConst) )
        if cap[e] > cap_max and cap[e] < 99:
          cap_max = cap[e]                  
    eLabel[e] = "{0:.1e}".format(cap[e])

# normalize the capacities to overcome difficulties with numerical precision 
for e in g.edges():
  if cap[e]< 99:
    cap[e] = cap[e]/cap_max                 # normalizes the capacities.

# solving the max flow problem
full_source = 0
full_target = 1
src, tgt = g.vertex(full_source), g.vertex(full_target)               # set the source and target nodes
res = gt.edmonds_karp_max_flow(g, src, tgt, cap)  
part = gt.min_st_cut(g, src, cap, res)
mc = sum([cap[e] - res[e] for e in g.edges() if part[e.source()] != part[e.target()]])
res.a = cap.a - res.a  # the actual flow
max_flow = sum(res[e] for e in tgt.in_edges())    # the maximum flow is the sum of the individual paths
flowLabel = g.new_edge_property("string") # initialize label for the edges in the graph

if abs(max_flow - 0.0)<eps:
  print ""
  print "Graph not connected. No percolation possible."
  print ""
  sys.exit()

q= 0
for v in g.vertices():  
  q=q+1

print q



try: # plot GEOS flux and pressure values on edges and nodes
  csv = np.genfromtxt (V_PressureFile, delimiter=";")
  P = csv[1:,1]
  vert_P = csv[1:,2]
  ok = 1
  maxP = max(P) 
  minP = min(P)
  Prel = (P-minP)/(maxP-minP)
  Prel[Prel<0.0]=0.0
  v_P              = g.new_vertex_property("double") 
  v_Prel           = g.new_vertex_property("double") 
  for i in xrange(len(P)):
    v = g.vertex(vert_P[i])
    v_P[v]    = P[i]
    v_Prel[v] = Prel[i] 

  if version == "seg" or version == "seg002":
    csvFlux  = np.genfromtxt (E_Fluxfile, delimiter=";")
    edgeFlux = csvFlux[1:,1]
    eSource  = csvFlux[1:,2]
    eTarget  = csvFlux[1:,3]
    e_max = len(edgeFlux)
    print "Flux extracted from file."

  print "Pressure extracted from file."
  print ""
 
except Exception as e:
  print "No pressure file found."
  ok = 0


v_pressure_color = g.new_vertex_property("vector<float>")
if ok == 1:
  for j,v in enumerate(g.vertices()):
    if j >1:     
      v_pressure_color[v] = [v_Prel[v],0.0, 1.0- v_Prel[v], 1.0]
      vLabel[v] = "{0:1.1f}".format(v_P[v]/1.0e6)
      print v_P[v]
    elif j == 0:
      v_pressure_color[v] = [1.0, 0.0, 0.0, 1.0]
      
    else:
      v_pressure_color[v] = [0.0, 0.0, 1.0, 1.0]

else:
  v_pressure_color = [0.0, 1.0, 0.0, 1.0]

e_penwidth = g.new_edge_property("double")
e_color    = g.new_edge_property("vector<float>")
resCap     = cap.copy()
e_Qp       = g.new_edge_property("double")            # GEOS volumerate result using pressures on vertices and cubic law.
e_Qv       = g.new_edge_property("double")            # GEOS volumerate result using velocities over the fluxAreabox.

if ok == 1 and (version == "seg" or version == "seg002") :
  for i in xrange(e_max):
    e = g.edge(eSource[i], eTarget[i])
    e_Qv[e] = edgeFlux[i]
    if abs(width[e] - 0.02) < eps:
      e_Qv[e] = e_Qv[e]/(0.4/diff)    # interpolates the arbitrary value, because the full values of the two adjacent cells are found here.

if normalize==0:
  for i, e in enumerate(g.edges()):
    if ok == 1:
      src = int(e.source()) 
      tgt = int(e.target())
      dP  = v_P[tgt] - v_P[src]                     # pressure difference between the two nodes.
      
      if length[e] > eps:                               # only for those edges that have a non-zero length calculate the volumerate approximations.
        e_Qp[e] = cap[e]/length[e]*abs(dP)*cap_max

      else:
        e_Qp[e] = 0.0
        
    else:
      e_Qp[e] = 0.0
      
    opacity = (0.001*i)
    
    if res[e] < 99 and res[e]>eps:
      flowLabel[e] = "{0:2.1f}%".format(res[e]/cap[e]*100)  # assign flow value to each edge while undoing the normalization

      if version == "seg" or version == "seg002":
        eLabel[e]    = flowLabel[e] +  " / "  +  "{0:.1f}".format(cap_max*cap[e]/1.0e-13) + " / " + "{0:.1f}".format(e_Qv[e]/1.0e-8)

      else:
        eLabel[e]    = flowLabel[e] +  " / "  +  "{0:.1f}".format(cap_max*cap[e]/1.0e-13)

      e_color[e]   = [0.0,0.0,0.0,1.]
      resCap[e]    = res[e] 

    elif res[e]>99:
      flowLabel[e] = "{0:.1e}".format(res[e])         # used to check for incorrect assumptions. If something is wrong, the graph will have 100.0 as the maximum flow on these edges.
      eLabel[e]    = flowLabel[e] +  " / "  +  eLabel[e]
      e_color[e] = [0.0,0.0,0.0,1.0]
      resCap[e]    = res[e] 

    else:
      eLabel[e] = ""
      resCap[e]    = cap[e]

elif normalize==1:
  for e in g.edges():
    src = e.source()
    tgt = e.target()
    flowLabel[e] = str(res[e])         # assign flow value to each edge while keeping the normalization
    if src>tgt:
      eLabel[e] = ""
    eLabel[e]    = eLabel[e] + " / " + flowLabel[e]

q= 0
for v in g.vertices():
  q=q+1

print q

##############################################################3
### Parameters used to draw simple algorithms.
##############################################################3
gt.graph_draw(g , pos               = pos,  # This draws the capacity graph using the positions and labels described above. The width of an edge is proportional to the capacity of that edge.
                  vertex_text       = vLabel,
                  vertex_font_size  = 20,
                  vertex_size       = 2, 
                  vertex_fill_color = v_pressure_color,
                  vertex_shape      = "square",
                  edge_text         = eLabel,
                  edge_font_size    = 14,
                  edge_start_marker = "none",
                  edge_end_marker   = "none",
                  edge_color        = e_color,
                  edge_text_distance = 26,
                  edge_pen_width    = gt.prop_to_size(resCap, mi=1.0, ma=10.0, log=False, power=0.1),
                  fmt               = "png",
                  fit_view          = True,
                  output            = outputPicGraph,
                  bg_color          = [1,1,1,1], # make the background white.
                  output_size       = (2500,1000) )
##############################################################3
##############################################################3

gt.graph_draw(g, pos               = pos, # This draws the flow graph using the positions and labels described above.
                 vertex_text       = vLabel, 
                 edge_pen_width    = gt.prop_to_size(res, mi=1 , ma=10, log=False , power=0.5),
                 edge_text         = flowLabel,
                 vertex_font_size  = 10,
                 vertex_size       = 1,
                 fit_view          = True,
                 vertex_fill_color = part,
                 edge_font_size    = 10,
                 fmt               = "png",
                 output            = outputpicMCFlow,
                 bg_color          = [1,1,1,1], # make the background white.
                 output_size       = (3000,3000) ) 


print "Max flow result:"
if normalize==0:
  print max_flow*cap_max # undoing the normalizaiton.
elif normalize==1:
  print max_flow         # keeping the normalization.

print "remove edges without flow"
noFlow = gt.find_edge_range(g, res, [0.0, eps])
for s in noFlow:
  g.remove_edge(s)




############################################################################################
### Get path information
############################################################################################
if ok == 0:
  g = gt.GraphView(g,vfilt=lambda x: x.in_degree()+x.out_degree()>0) # removes any remaining no edge vertices.

gCopy = g.copy() # This copy of g will be used to remove the flow in individual paths.

posCopy       = gCopy.vertex_properties["pos"]      # extracts node positions
resCopy       = gCopy.new_edge_property("double")      # capacity for flow in each edge
e_colorCopy   = gCopy.new_edge_property("vector<float>")

counter = 0
for v in gCopy.vertices():
  posCopy[v] = pos[v]
  counter+=1

print "amount of vertices", counter
counter = 0
for e in gCopy.edges():  
  resCopy[e] = res[e]
  e_colorCopy[e] = e_color[e]
  counter+=1

print "amount of edges", counter

if not os.path.exists(workingDir + "Paths"):
          os.makedirs(workingDir + "Paths")

def pathFinder(gCopy, src, tgt, resCopy, cap):
  pathEdges    = []
  pathLength   = []
  pathWidth    = []
  pathFlow     = []
  pathAperture = []
  pathCap      = []
  pathFactors  = []
  pathLdiag    = []
  pathGrav     = []
  numPaths     = 0
  resFull      = resCopy.copy()
  
  
  
  for i in xrange(0,1000): # perhaps do a while loop instead. This should be safer...
    pLength    = []  # store lengths of this path
    pWidth     = []  # store widths of this path
    pAperture  = []  # store apertures of this path
    pCap       = []  # store capacities of this path
    pFactor    = []  # store flow factor of this path (percentage of the edge that this path uses).
    pRes       = []  # store the flow results for this path
    

    vShort, eShort = tp.shortest_path(gCopy,src,tgt) # find shortest path
    
    if len(eShort) == 0: # no more path found
      print "no more paths"
      print "number of paths found: ", numPaths
      print i
      return pathEdges, pathLength, pathWidth, pathFlow, pathAperture, pathCap, numPaths, pathFactors, pathLdiag, pathGrav  # return values
    
    else:
      minFlow = 10000
      
      for e in eShort:
        e_colorCopy[e] = [1.0,0.0,0.0,1.0]
        pRes.append(resFull[e])     
        if abs(width[e]) > eps: 
          pWidth.append(width[e])
          pLength.append(length[e])
          pCap.append(cap[e])


        if resCopy[e]<minFlow:
          minFlow = resCopy[e]
      
      pFactor = minFlow/np.array(pRes)
       
      plot_path = 1
      if plot_path == 1:
          outputpicPath = workingDir + "Paths/" + inputFile.split("/")[3].split(".")[0] + "_Path_{0}.png".format(i+1)
          
          gt.graph_draw(gCopy , pos     = posCopy,  # This draws the capacity graph using the positions and labels described above. The width of an edge is proportional to the capacity of that edge.
                      vertex_font_size  = 0,
                      vertex_size       = 4, 
                      vertex_fill_color = v_pressure_color,
                      edge_font_size    = 0,
                      edge_color        = e_colorCopy,
                      edge_pen_width    = gt.prop_to_size(resCopy, mi=1.0, ma=2.0, log=False, power=0.1),
                      fmt               = "png",
                      output            = outputpicPath,
                      bg_color          = [1,1,1,1], # make the background white.
                      output_size       = (1000,1000) )
          sys.exit()
      
      for e in eShort:
        resCopy[e] = resCopy[e] - minFlow
        e_colorCopy[e] = [0.0,0.0,0.0,1.0]
      v_in   = vShort[1]
      v_out  = vShort[-2]
      pLdiag = np.sqrt( (cent[v_out][0] - cent[v_in][0])**2 + (cent[v_out][1] - cent[v_in][1])**2 + (cent[v_out][2] - cent[v_in][2])**2) # radial distance between inflow and outflow 
      pGrav  = 1.0e6 - 9.81*1000.0*(cent[v_out][2] - cent[v_in][2])

      noFlow = gt.find_edge_range(gCopy, resCopy, [0.0, eps])
      for s in noFlow:
        gCopy.remove_edge(s)
      if len(pWidth) >0:
        pathEdges.append([])
        pathEdges[numPaths]     = eShort
        pathLength.append([])
        pathLength[numPaths]    = pLength
        pathWidth.append([])
        pathWidth[numPaths]     = pWidth
        pathFlow.append([])
        pathFlow[numPaths]      = minFlow
        pathAperture.append([])
        pathAperture[numPaths]  = pAperture
        pathCap.append([])
        pathCap[numPaths]       = pCap
        pathFactors.append([])
        pathFactors[numPaths]   = pFactor
        pathLdiag.append(pLdiag)
        pathGrav.append(pGrav)
        numPaths += 1


src, tgt = g.vertex(full_source), g.vertex(full_target)               # set the source and target nodes
pathEdges, pathLength, pathWidth, pathFlow, pathAperture, pathCap, numPaths, pFactor, pathLdiag, pathGrav = pathFinder(gCopy, src, tgt, resCopy, cap)

pathLdiag = np.array(pathLdiag)



wMin  = []
wMax  = []
wHarm = []
wSuperHarm = []
Ltot  = []
Lmin  = []
Lmax  = []
capHarm = []
for n in xrange(numPaths):
  wMin.append(np.min(pathWidth[n]))
  wMax.append(np.max(pathWidth[n]))
  Ltot.append(sum(pathLength[n]))
  Lmin.append(np.min(pathLength[n]))
  Lmax.append(np.max(pathLength[n]))

print Ltot

for i, path in enumerate(pathLength):
  L = Ltot[i]
  wWeighted       = 0.0
  wSuperWeighted  = 0.0
  capWeighted     = 0.0
  eFactor = pFactor[i]
  
  for j,Li in enumerate(path):
      wWeighted = wWeighted + Li/L*1.0/pathWidth[i][j]
      wSuperWeighted = wSuperWeighted + Li/L*1.0/(pathWidth[i][j]*eFactor[j])
      capWeighted = capWeighted + Li/L*1.0/pathCap[i][j]
  
  
  wHarm.append(wWeighted**(-1))
  wSuperHarm.append(wSuperWeighted**(-1))
  capHarm.append(capWeighted**(-1))
 
wMin = np.array(wMin)


wHarm       = np.array(wHarm)
wSuperHarm  = np.array(wSuperHarm)
capHarm     = np.array(capHarm)

Ltot = np.array(Ltot)
a=1.0e-5
Kmin  = wMin*a**3/(12*mu*Ltot)
Kharm = wHarm*a**3/(12*mu*Ltot)
KSuperharm = wSuperHarm*a**3/(12*mu*Ltot)
KSuperharm_short = wSuperHarm*a**3/(12*mu*10.6)
KSuperharm_short2 = wSuperHarm*a**3/(12*mu*10.6)*Ltot/10.6
KSuperharm_diagonal = wSuperHarm*a**3/(12*mu*np.sqrt(10.6*10.6 + 6.0*6.0))
KSuperharm_diagonal2 = wSuperHarm*a**3/(12*mu*pathLdiag)
KSuperharm_diagonal3 = wSuperHarm*a**3/(12*mu*np.sqrt(10.6*10.6 + 6.0*6.0+ 6.0*6.0))*Ltot/np.sqrt(10.6*10.6 + 6.0*6.0 + 6.0*6.0)
KSuperharm_diagonal3 = wSuperHarm*a**3/(12*mu*np.sqrt(10.6*10.6 + 8.0*8.0+ 6.0*6.0))#*Ltot/np.sqrt(10.6*10.6 + 9.0*9.0)
KSuperharm_short3 = wSuperHarm*a**3/(12*mu*Ltot)*10.6

Qgrav = KSuperharm_diagonal2*np.array(pathGrav)
KminTot = sum(Kmin)
print "KminTot"
print KminTot

KharmTot = sum(Kharm)
print "KharmTot"
print KharmTot
print ""

KSuperharmTot = sum(KSuperharm)
print "KSuperharmTot"
print KSuperharmTot
print ""

KSuperharm_shortTot = sum(KSuperharm_short)
print "KSuperharm_shortTot"
print KSuperharm_shortTot
print ""

KSuperharm_shortTot2 = sum(KSuperharm_short2)
print "KSuperharm_shortTot2"
print KSuperharm_shortTot2
print ""


KSuperharm_shortTot3 = sum(KSuperharm_short3)
print "KSuperharm_shortTot3"
print KSuperharm_shortTot3
print ""


KSuperharm_diagonalTot = sum(KSuperharm_diagonal)
print "KSuperharm_diagonalTot"
print KSuperharm_diagonalTot
print ""

KSuperharm_diagonalTot2 = sum(KSuperharm_diagonal2)
print "KSuperharm_diagonalTot2"
print KSuperharm_diagonalTot2
print ""

KSuperharm_diagonalTot3 = sum(KSuperharm_diagonal3)
print "KSuperharm_diagonalTot3"
print KSuperharm_diagonalTot3
print ""

QgravTot = sum(Qgrav)
print "QgravTot"
print QgravTot
print ""

print "number of paths found: ", numPaths
print ""

total_flow = np.sum(pathFlow)
print total_flow*cap_max

if abs(max_flow - sum(pathFlow)) <= 1.0e-9:
      print "max flow result equal to path flow result"
else:
      print "flow amount error!!"

print cap_max
