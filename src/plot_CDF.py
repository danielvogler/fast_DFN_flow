import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#rc('text', usetex=True)

simulation = 0

small = 0
large = 1000
simName = "bm110_f250/"
workingDir    = "../benchmark/" + simName 
Qdscsv = np.genfromtxt (workingDir + "/Output/" + "results.csv", delimiter=",")
rhoFrac  = np.genfromtxt (workingDir + "rhoFrac.txt", delimiter=",")
Qhspmcsv = np.genfromtxt (workingDir + "paper_CDF_HSPM_investigation_results_Quartz.txt", delimiter=",")
Qispmcsv = np.genfromtxt (workingDir + "paper_CDF_ISPM_investigation_results_Quartz.txt", delimiter=",")

Qds   = []
Qhspm = []
Qispm = []

print( Qdscsv )

for simulation in range(small+1, large+1):
	Qds  .append(Qdscsv[simulation,1])
	Qhspm.append(Qhspmcsv[simulation,1])
	Qispm.append(Qispmcsv[simulation,1])

rho = np.sort(rhoFrac)




Hspmcounter = 0
Ispmcounter = 0
for i in range(len(Qhspm)):
	if Qhspm[i] == 0:
		Hspmcounter += 1
	if Qispm[i] == 0:
		Ispmcounter += 1
print( "# non-percolating DFN according to HSPM", Hspmcounter )
print( "# non-percolating DFN according to ISPM", Ispmcounter )
Qds   = np.array(Qds)  *1.0e6
Qhspm = np.array(Qhspm)*1.0e6
Qispm = np.array(Qispm)*1.0e6

Q1 = np.sort(Qds)
Q2 = np.sort(Qhspm)
Q3 = np.sort(Qispm)



N = len(Q2)
assert len(Q2) == N
assert len(Q3) == N
F1 = np.array(range(N))/float(N)
F2 = np.array(range(N))/float(N)
F3 = np.array(range(N))/float(N)
assert len(F1) == N
assert len(F2) == N
assert len(F3) == N



# Calculate D-values for the mesh refinement
D = 0.0021001239*rhoFrac - 0.0821045236
dOverR = 0.2/1.45
refinementFactor = 1+D*dOverR
QdsRefined = Qds/refinementFactor
Q1refined  = np.sort(QdsRefined)




Qmin = np.min([np.min(Q1refined),np.min(Q2),np.min(Q3)])

Q1MPS = Q1refined - Qmin*np.ones(len(Q1))
Q2MPS = Q2        - Qmin*np.ones(len(Q2))
Q3MPS = Q3        - Qmin*np.ones(len(Q3))

QMPSmax = np.max([np.max(Q1MPS),np.max(Q2MPS),np.max(Q3MPS)])
Q1MPS = Q1MPS/QMPSmax
Q2MPS = Q2MPS/QMPSmax
Q3MPS = Q3MPS/QMPSmax



colors = ['darkgreen', 'yellowgreen', 'springgreen', 'thistle', 'teal', 'seagreen', 'palegreen', 'olivedrab', 'olive', 'mediumspringgreen', 'mediumseagreen', 'limegreen', 'lightseagreen', 'lightgreen', 'greenyellow', 'green', 'forestgreen', 'darkseagreen', 'darkolivegreen', 'darkgreen']


cDS   = 'k'          # DS   line color
cHSPM = '#99ff66'  # HSPM line color 
cISPM = 'darkorange' # ISPM line color
xyFontsize   = 21  	 # fontsize for axis   title
tFontsize    = 16    # fontsize for figure title
tickFontsize = 20    # fontsize for axis labels
lFontsize    = 20	 # fontsize for the legend

plt.rcParams['xtick.major.pad']='10'
nFig = 1
plt.figure(nFig)
plt.plot(Q1MPS,F1, '-', color=cDS, linewidth=2)
plt.plot(Q2MPS,F2, '-', color=cHSPM, linewidth=2)
plt.plot(Q3MPS,F3, '-', color=cISPM, linewidth=2)
plt.xlabel("Q (normalized)", fontsize=xyFontsize)
plt.ylabel("CDF"           , fontsize=xyFontsize)
plt.legend(["Q$_{DS}$", "Q$_{HSPM}$", "Q$_{ISPM}$"], loc="lower right", fontsize = lFontsize)
plt.tick_params(axis='both', which='major', labelsize=tickFontsize)
plt.xticks(np.arange(0.0, 1.5, 0.5))
plt.yticks(np.arange(0.0, 1.5, 0.5))
plt.tight_layout()
nFig += 1



binsize = int(np.sqrt(large).round())
print("binsize", binsize)
binsize = 16


binwidth = 1.0/binsize
plt.figure(nFig)
plt.hist(Q1MPS, bins=np.arange(0.0, 1.0 + binwidth, binwidth), histtype = 'step', color=cDS  , label='Q$_{DS}$', linewidth=2)
plt.hist(Q2MPS, bins=np.arange(0.0, 1.0 + binwidth, binwidth), histtype = 'step', color=cHSPM, alpha=1.0, label="Q$_{HSPM}$", linewidth=2)
plt.hist(Q3MPS, bins=np.arange(0.0, 1.0 + binwidth, binwidth), histtype = 'step', color=cISPM, alpha=1.0, label='Q$_{ISPM}$', linewidth=2)
plt.xlabel("Q (normalized)"             , fontsize=xyFontsize)
plt.ylabel("Realizations"               , fontsize=xyFontsize)
plt.legend(fontsize = lFontsize-2, loc='upper right')
plt.tick_params(axis='both', which='major', labelsize=tickFontsize)
plt.xticks(np.arange(0.0, 1.5, 0.5))
plt.yticks(np.arange(0.0, 240, 80))
plt.tight_layout()
nFig += 1


plt.show()



