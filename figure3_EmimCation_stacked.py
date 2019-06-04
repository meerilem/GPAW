import matplotlib.pyplot as plt
import numpy as np

numberlist = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
# 
def readTheFile(filename):
	m, atomname = [], []
	notatom=True
	with open(filename) as f:
		count = 0
		for line in f:
			count +=1
			if count > 1 and count%2 == 1:
				m.append(float(line[:-1]))
			elif count > 1 and count%2 == 0:
			    tmp = line.split()
			    atomname.append(tmp[0])
	return m, atomname
					
def gaussian(x, m, s):
  return np.exp(-np.power((x - m), 2.) / (2 * np.power(s, 2.)))

def getListOfAtoms(atom, m, atomname):
	m_atom, atomname_atom = [], []
	for i in range(len(atomname)):
	    if len(atom) == 1:
	        if atomname[i][0] == atom and atomname[i][1] in numberlist:
	            m_atom.append(m[i])
	            atomname_atom.append(atomname[i])
	    else:
	        if atomname[i][0:2] == atom:
	            m_atom.append(m[i])
	            atomname_atom.append(atomname[i])
	
	return m_atom, atomname_atom

def plotAtomSpectra(m, atomname, s):
	t = []
	w = np.linspace(min(m) - 1, max(m) + 1, 201)
	for i in range(len(m)):
	    y_pos = 1.1
	    for j in range(i):
	        if abs(m[i] - m[j]) < 0.3:
	            y_pos += 0.1	        
	for j in range(len(w)):
	  a = 0
	  for i in range(len(m)):
	    a += gaussian(w[j], m[i], s)
	  t.append(a)
	return w, t
	
def plotFinalSpectra(namelist, nr, lineW, dashes, ncolors, plot_list):
	w, t = [], []
	count = 0
	for name in namelist:
		m, atomname = [], []
		m, atomname = readTheFile('./ready/' + name + '/' + name + filename + '.out')
		m_tmp, atomname = getListOfAtoms(atom, m, atomname)
		m2 = [x+aliphatic-m_tmp[0] for x in m_tmp]
		w_tmp, t_tmp = plotAtomSpectra(m2, atomname , s)
		w.append(w_tmp)
		new_list = [x+nr[count]*4 for x in t_tmp]
		t.append(new_list)
		count += 1

	for i in range(len(namelist)):
		plot_list += plt.plot(w[i], t[i], linewidth=lineW, c=ncolors[i], label=namelist[i], dashes=dashes)

########################################################
ncolors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']
black=['k'] * 8
grey=['#707070'] * 8

s = 0.5
atom = 'C'

ax = plt.figure(figsize=(8.3/2.54,8.3/2.54)) # 8.3 cm to inches
butane = 289.65
aliphatic = 285.0

plot_list = []

namelist1 = ['EMImBCN4', 'EMImTFSI', 'EMImFSI','EMImPF6', 'EMImBF4', 'EMImCl', 'EMImBr', 'EMImI']
numbers = [0, 1, 2, 3, 4, 5, 6, 7]
filename = ""
plotFinalSpectra(namelist1, numbers, 3.0, [1,0], ncolors,plot_list)

	
namelist2 = ['EMImTFSI','EMImPF6', 'EMImBF4', 'EMImCl', 'EMImBr', 'EMImI']
nr = [1, 3, 4, 5, 6, 7]
filename = "_Villar"
plotFinalSpectra(namelist2, nr, 1.0, [4,1,4,1], black,plot_list)

namelist2 = ['EMImBCN4']
nr = [0]
filename = "_Kruusma"
plotFinalSpectra(namelist2, nr, 1.0, [4,1,1,1], black,plot_list)

	
namelist2 = ['EMImBF4']
nr = [4]
filename ='_Tonisoo'
plotFinalSpectra(namelist2, nr, 1.0, [4,1,1,1], grey,plot_list)

namelist2 = ['EMImBCN4', 'EMImTFSI', 'EMImBF4']
nr = [0, 1, 4]
filename = '_Kotz'
plotFinalSpectra(namelist2, nr, 1.0, [1,1], black,plot_list)

namelist2 = ['EMImTFSI']
nr = [1]
filename = '_Hammer'
plotFinalSpectra(namelist2, nr, 1.0, [4,1,4,1,1,1], grey,plot_list)

namelist2 = ['EMImTFSI']
nr = [1]
filename = '_Reinmoller'
plotFinalSpectra(namelist2, nr, 1.0, [4,1,4,1,1,1,1,1], black,plot_list)

"""
from matplotlib.lines import Line2D
legend_elements = []
for i in range(len(namelist1)):
	legend_elements.append(Line2D([0], [0], ls='-', color=ncolors[7-i], label=namelist1[7-i], markersize=10))
#legend_elements.append(Line2D([0], [0], ls='-', color='k', label='calc', markersize=10))
legend_elements.append(Line2D([0], [0], ls='--', color='grey', label='Villar', markersize=10))
legend_elements.append(Line2D([0], [0], ls='-.', color='grey', label='Tonisoo', markersize=10))
legend_elements.append(Line2D([0], [0], ls=':', color='k', label='Schmitz', markersize=10))
legend_elements.append(Line2D([0], [0], ls='-', color='grey', label='Hammer', markersize=10))
legend_elements.append(Line2D([0], [0], ls='-.', color='k', label='Reinmoller', markersize=10))
legend_elements.append(Line2D([0], [0], ls='--', color='k', label='Kruusma', markersize=10))

plt.legend(handles=legend_elements, loc='best', prop={'size': 10})
"""

# * Loop to create multiple legends

from matplotlib.legend import Legend

# * Change legend labels and y-positions here

leg_names = [ r"$name_{1}$",r"$name_{2}$",
							r"$name_{3}$", r"$name_{4}$",
							r"$name_{5}$", r"$name_{6}$",
							r"$name_{7}$", r"$name_{8}$"]

leg_ypos = [0.22, 0.32, 0.405, 0.505,
						0.593, 0.68, 0.77, 0.875 ]

for i in range(8):
	leg = Legend(ax,[plot_list[i]],[leg_names[i]],loc=(0.73,leg_ypos[i]),frameon=False,fontsize=7,)
	ax.add_artist(leg)

# leg_1 = Legend(ax,plot_list[:2],[r'$G_{\mu \nu}$',r'$\Lambda g_{\mu \nu}$'],loc=(0.7,0.1),frameon=False)
# ax.add_artist(leg_1)

# leg_2 = Legend(ax,plot_list[2],['name3','name4'],loc=loc=(0.7,0.1),frameon=False)
# ax.add_artist(leg_2)

# leg_3 = Legend(ax,plot_list[4:6],['name5','name6'],loc='lower right',frameon=False)
# ax.add_artist(leg_3)

plt.yticks([])
plt.xticks([284,285,286,287,288])

#plt.title(r"%s %s1s XPS spectra" % ('EMIm cation', atom))
plt.ylabel("intensity")
plt.xlabel("binding energy / eV")
#plt.tick_params(labelsize=16)
plt.xlim(283.8,289.5)
plt.tight_layout()
plt.savefig('figure3_EmimCation_stacked.png', format="png", dpi=300, bbox_inches='tight')
plt.savefig('figure3_EmimCation_stacked.svg', format="svg", dpi=2000, bbox_inches='tight')
plt.savefig('figure3_EmimCation_stacked.eps', format="eps", dpi=2000, bbox_inches='tight')




