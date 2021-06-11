import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


simulation = 0
csv = np.genfromtxt ("100_DFN.csv", delimiter=",")

rhoFrac = []
Qds   	= []
Qhspm 	= []
Qispm 	= []
CPUds	= []
CPUhspm = []
CPUispm = []
edHSPM 	= []
edISPM  = []
Vhspm 	= []
Vispm 	= []


for i in range(1,101):
	rhoFrac.append(csv[i,0])
	Qds    .append(csv[i,1])
	Qhspm  .append(csv[i,3])
	Qispm  .append(csv[i,8])
	CPUds  .append(csv[i,2])
	CPUhspm.append(csv[i,5])
	CPUispm.append(csv[i,10])
	edHSPM .append(csv[i,6])
	edISPM .append(csv[i,11])
	Vhspm  .append(csv[i,7])
	Vispm  .append(csv[i,12])


rhoFrac = np.array(rhoFrac)
Qds 	= np.array(Qds 	)
Qhspm 	= np.array(Qhspm)
Qispm 	= np.array(Qispm)

rho 	= np.sort(rhoFrac)


Ehspm = (Qhspm/Qds - 1)*100
Eispm = (Qispm/Qds - 1)*100



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














cDS   = 'k'          # DS   line color
cHSPM = 'darkorange' # HSPM line color
cISPM = 'limegreen'  # ISPM line color
xyFontsize   = 21  	 # fontsize for axis   title
tFontsize    = 16    # fontsize for figure title
tickFontsize = 20    # fontsize for axis labels
lFontsize    = 20	 # fontsize for the legend

plt.rcParams['xtick.major.pad']='10'
nFig = 1







zds   = np.polyfit(rhoFrac, Qds  ,1)
zhspm = np.polyfit(rhoFrac, Qhspm,1)
zispm = np.polyfit(rhoFrac, Qispm,1)
pds   = np.poly1d(zds)
phspm = np.poly1d(zhspm)
pispm = np.poly1d(zispm)

plt.figure(nFig)
plt.plot(rhoFrac, pds(rhoFrac)  , cDS  , alpha=1.0)
plt.plot(rhoFrac, phspm(rhoFrac), cHSPM, alpha=1.0)
plt.plot(rhoFrac, Qds  , cDS  , linestyle = 'None', marker='s', alpha=1.0, fillstyle='full', label='Q$_{DS}$'  , markersize=10)
plt.plot(rhoFrac, Qhspm, cHSPM, linestyle = 'None', marker='^', alpha=1.0, fillstyle='full', label='Q$_{HSPM}$', markersize=10)
plt.xlabel("Number of Fractures" , fontsize=xyFontsize)
plt.ylabel("Q [mL/s]"            , fontsize=xyFontsize)
plt.tick_params(axis='both', which='major', labelsize=tickFontsize)
plt.xticks(np.arange(140, 360, 40))
plt.yticks(np.arange(0.1, 1.3, 0.4))
plt.legend(loc="upper left", fontsize = lFontsize-2)
plt.tight_layout()
nFig += 1


plt.figure(nFig)
plt.plot(rhoFrac, pds(rhoFrac)  , cDS  , alpha=1.0)
plt.plot(rhoFrac, pispm(rhoFrac), cISPM, alpha=1.0)
plt.plot(rhoFrac, Qds  , cDS  , linestyle = 'None', marker='s', alpha=1.0, fillstyle='full', label='Q$_{DS}$'  , markersize=10)
plt.plot(rhoFrac, Qispm, cISPM, linestyle = 'None', marker='o', alpha=1.0, fillstyle='full', label='Q$_{ISPM}$', markersize=10)
plt.xlabel("Number of Fractures" , fontsize=xyFontsize)
plt.ylabel("Q [mL/s]"            , fontsize=xyFontsize)
plt.tick_params(axis='both', which='major', labelsize=tickFontsize)
plt.xticks(np.arange(140, 360, 40))
plt.yticks(np.arange(0.1, 1.3, 0.4))
plt.legend(loc="upper left",fontsize = lFontsize-2,)
plt.tight_layout()
nFig += 1





plt.figure(nFig)
plt.plot(rhoFrac, Ehspm, cHSPM, linestyle = 'None', marker='^', alpha=1.0, fillstyle='full', label='Q$_{HSPM}$', markersize=10)
plt.xlabel("Number of Fractures" , fontsize=xyFontsize)
plt.ylabel("${\epsilon}$ [%]"            , fontsize=xyFontsize)
plt.tick_params(axis='both', which='major', labelsize=tickFontsize)
plt.xticks(np.arange(140, 360, 40))
plt.yticks(np.arange(-50, 30, 10))
plt.legend(loc="lower center", fontsize = lFontsize-2)
plt.tight_layout()
nFig += 1


plt.figure(nFig)
plt.plot(rhoFrac, Eispm, cISPM, linestyle = 'None', marker='o', alpha=1.0, fillstyle='full', label='Q$_{ISPM}$', markersize=10)
plt.xlabel("Number of Fractures" , fontsize=xyFontsize)
plt.ylabel("${\epsilon}$ [%]"            , fontsize=xyFontsize)
plt.tick_params(axis='both', which='major', labelsize=tickFontsize)
plt.xticks(np.arange(140, 360, 40))
plt.yticks(np.arange(-50, 30, 10))
plt.legend(loc="lower center", fontsize = lFontsize-2)
plt.tight_layout()
nFig += 1


zds   = np.polyfit(rhoFrac, CPUds  ,1)
zhspm = np.polyfit(rhoFrac, CPUhspm,1)
zispm = np.polyfit(rhoFrac, CPUispm,1)
pds   = np.poly1d(zds)
phspm = np.poly1d(zhspm)
pispm = np.poly1d(zispm)

plt.figure(nFig)
plt.semilogy(rhoFrac, pds(rhoFrac)  , cDS  , alpha=1.0)
plt.semilogy(rhoFrac, phspm(rhoFrac), cHSPM, alpha=1.0)
plt.semilogy(rhoFrac, pispm(rhoFrac), cISPM, alpha=1.0)
plt.semilogy(rhoFrac, CPUds  , cDS  , linestyle = 'None', marker='s', alpha=1.0, fillstyle='full', label='Q$_{DS}$'  , markersize=10)
plt.semilogy(rhoFrac, CPUhspm, cHSPM, linestyle = 'None', marker='^', alpha=1.0, fillstyle='full', label='Q$_{HSPM}$', markersize=10)
plt.semilogy(rhoFrac, CPUispm, cISPM, linestyle = 'None', marker='o', alpha=1.0, fillstyle='full', label='Q$_{ISPM}$', markersize=10)
plt.xlabel("Number of Fractures" , fontsize=xyFontsize)
plt.ylabel("Computational cost [s]"            , fontsize=xyFontsize)
plt.tick_params(axis='both', which='major', labelsize=tickFontsize)
plt.xticks(np.arange(140, 360, 40))
plt.legend(loc="lower center", fontsize = lFontsize-2)
plt.tight_layout()
nFig += 1

plt.figure(nFig)
plt.plot(rhoFrac, phspm(rhoFrac), cHSPM, alpha=1.0)
plt.plot(rhoFrac, pispm(rhoFrac), cISPM, alpha=1.0)
plt.plot(rhoFrac, CPUhspm, cHSPM, linestyle = 'None', marker='^', alpha=1.0, fillstyle='full', label='Q$_{HSPM}$', markersize=10)
plt.plot(rhoFrac, CPUispm, cISPM, linestyle = 'None', marker='o', alpha=1.0, fillstyle='full', label='Q$_{ISPM}$', markersize=10)
plt.xlabel("Number of Fractures" , fontsize=xyFontsize)
plt.ylabel("Computational cost [s]"            , fontsize=xyFontsize)
plt.tick_params(axis='both', which='major', labelsize=tickFontsize)
plt.xticks(np.arange(140, 360, 40))
plt.legend(loc="lower center", fontsize = lFontsize-2)
plt.tight_layout()
nFig += 1


zhspm = np.polyfit(rhoFrac, edHSPM,1)
zispm = np.polyfit(rhoFrac, edISPM,1)
pds   = np.poly1d(zds)
phspm = np.poly1d(zhspm)
pispm = np.poly1d(zispm)

plt.figure(nFig)
plt.plot(rhoFrac, phspm(rhoFrac), cHSPM, alpha=1.0)
plt.plot(rhoFrac, pispm(rhoFrac), cISPM, alpha=1.0)
plt.plot(rhoFrac, edHSPM, cHSPM, linestyle = 'None', marker='^', alpha=1.0, fillstyle='full', label='Q$_{HSPM}$', markersize=10)
plt.plot(rhoFrac, edISPM, cISPM, linestyle = 'None', marker='o', alpha=1.0, fillstyle='full', label='Q$_{ISPM}$', markersize=10)
plt.xlabel("Number of Fractures" , fontsize=xyFontsize)
plt.ylabel("Number of Edges"            , fontsize=xyFontsize)
plt.tick_params(axis='both', which='major', labelsize=tickFontsize)
plt.xticks(np.arange(140, 360, 40))
plt.legend(loc="lower right", fontsize = lFontsize-2)
plt.tight_layout()
nFig += 1

zhspm = np.polyfit(rhoFrac, Vhspm,1)
zispm = np.polyfit(rhoFrac, Vispm,1)
pds   = np.poly1d(zds)
phspm = np.poly1d(zhspm)
pispm = np.poly1d(zispm)

plt.figure(nFig)
plt.plot(rhoFrac, phspm(rhoFrac), cHSPM, alpha=1.0)
plt.plot(rhoFrac, pispm(rhoFrac), cISPM, alpha=1.0)
plt.plot(rhoFrac, Vhspm, cHSPM, linestyle = 'None', marker='^', alpha=1.0, fillstyle='full', label='Q$_{HSPM}$', markersize=10)
plt.plot(rhoFrac, Vispm, cISPM, linestyle = 'None', marker='o', alpha=1.0, fillstyle='full', label='Q$_{ISPM}$', markersize=10)
plt.xlabel("Number of Fractures" , fontsize=xyFontsize)
plt.ylabel("Number of Vertices"            , fontsize=xyFontsize)
plt.tick_params(axis='both', which='major', labelsize=tickFontsize)
plt.xticks(np.arange(140, 360, 40))
plt.legend(loc="lower right", fontsize = lFontsize-2)
plt.tight_layout()
plt.show()
nFig += 1


 
 
