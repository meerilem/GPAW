import matplotlib.pyplot as plt
import numpy as np


numberlist = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
atomname, m = [], []
ncolors = ['#ff0000', '#0000ff', '#00ff00', '#00ffff', '#ffd700', '#ff00ff', '#ff8000', '#80ff00', '#00ff80', '#0080ff', '#7f00ff', '#ff007f', '#bc8f8f', '#c71585', '#ff0000', '#0000ff', '#00ff00', '#00ffff', '#ffd700', '#ff00ff', '#ff8000', '#80ff00', '#00ff80', '#0080ff', '#7f00ff', '#ff007f', '#bc8f8f', '#c71585']
ncolors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']


def ReadChargeFile(filename):
	with open(filename, 'r') as f:
	    array = []
	    lines = f.readlines()
	    for line in lines:
	        atom.append(line.split()[0])
	        DDEC.append(line.split()[1])
	        bader.append(line.split()[3])
	        catan.append(line.split()[4])
            
	return atom, DDEC, bader, catan

def Charges(DDEC, bader, catan):
    DDEC_cat, DDEC_an = 0, 0
    bader_cat, bader_an = 0, 0
    for i in range(len(bader)):
        if catan[i] == 'AN':
            DDEC_an += float(DDEC[i])
            bader_an += float(bader[i])
        elif catan[i] == 'CAT':
            DDEC_cat += float(DDEC[i])
            bader_cat += float(bader[i])
                
    return DDEC_cat, DDEC_an, bader_cat, bader_an
    
def SortCharges(atom, DDEC, bader, whatAtom):
	DDEC_atom, bader_atom = [], []
	for i in range(len(DDEC)):
		if atom[i] == whatAtom:
			DDEC_atom.append(float(DDEC[i]))
			bader_atom.append(float(bader[i]))
	return DDEC_atom, bader_atom
			
	
def ReadSpectraFile(filename):
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
					

def GetListOfAtoms(atom):
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
	



########################################################
cation = ['TEPA', 'BPy', 'Pyr14', 'BMIm', 'EMIm']
anion = ['TFSI', 'FSI', 'BF4', 'PF6', 'BCN4', 'Cl', 'Br', 'I']
marker = ["o", "v", "^", "p", "8", "*", "s", "d", "+"] 
#anion = ['TFSI', 'FSI', 'BF4', 'PF6', 'BCN4']
#anion = [anion[7]]
cation = [cation[4]]
X, Y = [], []
nr_of_atom = 4
whatAtom = 'C'

plt.figure(figsize=(3.75, 3.75))

count = -1
for cat in cation:
	for an in anion:
	    count += 1
	    atom, DDEC, bader, catan, DDEC_atom, bader_atom = [], [], [], [], [], []
	    atom, DDEC, bader, catan = ReadChargeFile('./ready/' + cat + an + '/' + cat + an + '_charges.out')
	    DDEC_cat, DDEC_an, bader_cat, bader_an = Charges(DDEC, bader, catan)
	    DDEC_atom, bader_atom = SortCharges(atom, DDEC, bader, whatAtom)
	    for i in range(len(bader_atom)):
	    	nr_of_atom = i
	    	Y.append([bader_atom[nr_of_atom]])
	    	atomname, m = [], []
	    	m, atomname = ReadSpectraFile('./ready/' + cat + an + '/' + cat + an + '.out')
	    	X.append([m[nr_of_atom]])
	    plt.scatter(bader_atom, m[:len(bader_atom)], c=ncolors[count], marker = marker[count], s=16, label=an)
	 
from sklearn.linear_model import LinearRegression
reg = LinearRegression().fit(X, Y) 
y = reg.predict(X)
plt.plot(y, X, 'k-', lw=1.0)
plt.legend()
plt.xlabel("DDEC6 charge / e",fontsize=9)
plt.ylabel("binding energy / eV",fontsize=9)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)

plt.tight_layout()

plt.savefig('figure5_ChargeVSSpectra_%s%s_DDEC6_allAtoms.png' % (cation[0],"all"), dpi=2000, bbox_inches='tight')






