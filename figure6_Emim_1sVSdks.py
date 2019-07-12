import matplotlib.pyplot as plt
import numpy as np

numberlist = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
atomname, m = [], []
ncolors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']
ncolors=['#CC3333', '#CC3333', '#CC3333', '#CC3333', '#CC3333', '#CC3333', '#CC3333', '#0055ff']

def ReadTheFile1(filename):
	# notatom=True
	m, atomname = [], []
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

def GetListOfAtoms(atom, m, atomname):
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

def PlotAtomSpectra(m, atomname, s):
	t = []
	w = np.linspace(min(m) - 1, max(m) + 1, 201)
	# plot individual signals
	for i in range(len(m)):
	    #plt.plot(w, gaussian(w, m[i], s), label=atomname[i], color=ncolors[i])
	    y_pos = 1.1
	    for j in range(i):
	        if abs(m[i] - m[j]) < 0.3:
	            y_pos += 0.1	        
	    #plt.text(m[i], y_pos, atomname[i], horizontalalignment='center', rotation=0, size=10, color=ncolors[i])

	# plot total signal
	for j in range(len(w)):
	  a = 0
	  for i in range(len(m)):
	    a += gaussian(w[j], m[i], s)
	  t.append(a)

	return w, t

########################################################
cation = ['TEPA', 'BPy', 'Pyr14', 'BMIm', 'EMIm']
anion = ['TFSI', 'FSI', 'PF6', 'BCN4', 'Cl', 'Br', 'I', 'BF4']
# marker = ["o", "v", "^", "p", "8", "*", "s", "d", "+"] 
marker = ["x", "x", "x", "x", "x", "x", "x", "D"] 
alpha = [0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.9] 
#anion = [anion[0]]
cat = cation[4]
namelist = []
#an = 'Cl'	
atom='C'
s=0.25
count = -1

ax = plt.figure(figsize=(8.3/2.54,8.3/2.54)) 

X, Y = [], []
for an in anion:
	atomname, m = [], []
	m, atomname = ReadTheFile1('./data/' + cat + an + '/' + cat + an + '.out')	
	m_tmp, atomname = GetListOfAtoms(atom, m, atomname)
	m = [x for x in m_tmp]
	atomname1s, m1s = [], []
	m1s, atomname1s = ReadTheFile1('./data/' + cat + an + '/' + cat + an + '_1s.out')	
	m_tmp1s, atomname1s = GetListOfAtoms(atom, m1s, atomname1s)
	m1s = [x1 for x1 in m_tmp1s]
	count += 1
	for i in range(len(m1s)):
	    nr_of_atom = i
	    Y.append([m[nr_of_atom]])
	    X.append([m1s[nr_of_atom]])
	if an == 'BF4':
		bf4 = plt.scatter(m, m1s, c=ncolors[count], marker = marker[count], s=16, label=an,alpha=alpha[count],edgecolors="#333366",linewidths=0.7)
	else:
		plt.scatter(m, m1s, c=ncolors[count], marker = marker[count], s=16, label=an,alpha=alpha[count],edgecolors="#333366",linewidths=0.7)
	 
from sklearn.linear_model import LinearRegression
reg = LinearRegression().fit(X, Y) 
xplot = np.linspace(start=268,stop=274,num=100).reshape(-1,1)
y = reg.predict(xplot)
plt.plot(y, xplot, 'k-', lw=1.0,linestyle='--')
params = {'legend.fontsize': 7}
plt.rcParams.update(params)


# from matplotlib.legend import Legend

# leg = Legend(ax,[bf4],[r"BF$_{4}$"],loc=(0.22,0.85),
# 	frameon=False,fancybox=True, fontsize=7,handlelength=0)
# ax.add_artist(leg)

plt.xlim(left=288.8, right=293.8)
plt.ylim(bottom=268.2, top=273.2)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.minorticks_on()
plt.xlabel(r"$\Delta$" + "KS binding energy / eV", fontsize=8)
plt.ylabel("1s binding energy / eV", fontsize=8)
plt.tight_layout()
plt.savefig('./figures_for_article/figure6_Spectra1sVSSpectradKs_%s%s.png' % (cation[0],"all"), dpi=300, bbox_inches='tight')
plt.savefig('./figures_for_article/figure6_Spectra1sVSSpectradKs_%s%s.svg' % (cation[0],"all"), dpi=600, bbox_inches='tight')
#plt.savefig('./figures_for_article/figure6_Spectra1sVSSpectradKs_%s%s.eps' % (cation[0],"all"), dpi=2000, bbox_inches='tight')
