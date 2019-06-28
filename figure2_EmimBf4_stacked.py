import matplotlib.pyplot as plt
import numpy as np

numberlist = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
ncolors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']

def readTheFile(filename):
	m, atomname = [], []
	# notatom=True
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

def plotAtomSpectra(m, atomname, s, y_position, dashes):
	atoms = ['C6', 'C5', 'C1', 'C2', 'C3', 'C4', 'X', 'Y']
	t = []
	w = np.linspace(min(m) - 3, max(m) + 3, 201)
	# plot individual signals
	for i in range(len(m)):
	    plt.plot(w, gaussian(w, m[i], s) + y_position, label=atomname[i], color=ncolors[i], dashes=dashes)
	    y_pos = 1.1
	    for j in range(i):
	        if abs(m[i] - m[j]) < 0.1:
	            y_pos += 0.3	        
	    plt.text(m[i], y_pos + y_position, atoms[i], horizontalalignment='center', rotation=0, size=7, color=ncolors[i])
	# plot total signal
	for j in range(len(w)):
	  a = y_position
	  for i in range(len(m)):
	    a += gaussian(w[j], m[i], s)
	  t.append(a)
	return w, t
	
def plotFinalSpectra(lineW, dashes, y_pos):
	w, t = [], []
	m, atomname = [], []
	m, atomname = readTheFile('./data/' + name + '/' + name + filename + '.out')
	m_tmp, atomname = getListOfAtoms(atom, m, atomname)
	m2 = [x+aliphatic-m_tmp[0] for x in m_tmp]
	w, t = plotAtomSpectra(m2, atomname , s, y_pos, dashes)
	plt.plot(w, t, linewidth=lineW, color='k', dashes=dashes)	
########################################################
s = 0.5
atom = 'C'

fig = plt.figure(figsize=(8.3/2.54,8.3/2.54))
#fig, ax = plt.subplots()

butane = 289.65
aliphatic = 285.0

linewidth = 1

name = 'EMImBF4'
atom = 'C'
filename = ""
plotFinalSpectra(linewidth, [1,0], 0)

filename='_Villar'
plotFinalSpectra(linewidth, [3,1,3,1], 2)

filename='_Tonisoo'
plotFinalSpectra(linewidth, [3,1,1,1], 4)

filename='_Kotz'
plotFinalSpectra(linewidth, [1,1], 6)

filename='_1s'
plotFinalSpectra(linewidth, [1,1,1,1,3,1], -2)

from matplotlib.lines import Line2D
legend_elements = []
legend_elements.append(Line2D([0], [0], dashes=[1,1], color='k', label=r"$Schmitz$", markersize=10))
legend_elements.append(Line2D([0], [0], dashes=[3,1,1,1], color='k', label=r"$Tonisoo$", markersize=10))
legend_elements.append(Line2D([0], [0], dashes=[3,1,3,1], color='k', label=r"$Garcia$", markersize=10))
legend_elements.append(Line2D([0], [0], dashes=[1,0], color='k', label=r"${\Delta}KS$", markersize=10))
legend_elements.append(Line2D([0], [0], dashes=[1,1,1,1,3,1], color='k', label=r"$1s$", markersize=10))
params = {'legend.fontsize': 7,
          'legend.handlelength': 2}
plt.rcParams.update(params)
plt.legend(handles=legend_elements, loc='upper right')
plt.xlim(283.5,289)
plt.ylabel("intensity / arb. units",fontsize=8)
plt.xlabel("binding energy / eV",fontsize=8)
plt.xticks(fontsize=7)
plt.yticks([])

from matplotlib.cbook import get_sample_data
path="EMIm.png"
im = plt.imread(path)

newax = fig.add_axes([0.12, 0.62, 0.3, 0.3], anchor='NE', zorder=1)
plt.imshow(im)
newax.imshow(im)
newax.axis('off')
#newax.axes.set_yticklabels([])

plt.xticks(fontsize=7)
plt.yticks([])
#plt.title(r"%s %s1s XPS spectra" % ('EMImBF4', atom))
plt.ylabel("intensity / arb. units",fontsize=8)
plt.xlabel("binding energy / eV",fontsize=8)
fig.tight_layout()
fig.savefig('./figures_for_article/figure2_EmimBf4_stacked.png', format="png", dpi=300, bbox_inches='tight')
fig.savefig('./figures_for_article/figure2_EmimBf4_stacked.svg', format="svg", dpi=2000, bbox_inches='tight')
fig.savefig('./figures_for_article/figure2_EmimBf4_stacked.eps', format="eps", dpi=2000, bbox_inches='tight')

