import matplotlib.pyplot as plt
import numpy as np


numberlist = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
atomname, m = [], []
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

plt.figure(figsize=(8.3/2.54,8.3/2.54)) # 8.3 cm to inches

count = -1
for cat in cation:
	for an in anion:
	    count += 1
	    atom, DDEC, bader, catan, DDEC_atom, bader_atom = [], [], [], [], [], []
	    atom, DDEC, bader, catan = ReadChargeFile('./data/' + cat + an + '/' + cat + an + '_charges.out')
	    DDEC_cat, DDEC_an, bader_cat, bader_an = Charges(DDEC, bader, catan)
	    DDEC_atom, bader_atom = SortCharges(atom, DDEC, bader, whatAtom)
	    atomname, m = [], []
	    m, atomname = ReadSpectraFile('./data/' + cat + an + '/' + cat + an + '.out')
	    for i in range(len(bader_atom)):
	    	nr_of_atom = i
	    	Y.append([bader_atom[nr_of_atom]])
	    	X.append([m[nr_of_atom]])
	    plt.scatter(bader_atom, m[:len(bader_atom)], c=ncolors[count], marker = marker[count], s=16, label=an)
	 
from sklearn.linear_model import LinearRegression
reg = LinearRegression().fit(X, Y) 
y = reg.predict(X)
plt.plot(y, X, 'k-', lw=1.0)
params = {'legend.fontsize': 7}
plt.rcParams.update(params)
plt.legend()
plt.xlabel("DDEC6 charge / e",fontsize=8)
plt.ylabel(r"$\Delta$" + "KS binding energy / eV",fontsize=8)
plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.tight_layout()
plt.savefig('./figures_for_article/figure5_ChargeVSSpectra_%s%s_DDEC6.png' % (cation[0],"all"), format="png", dpi=300, bbox_inches='tight')
plt.savefig('./figures_for_article/figure5_ChargeVSSpectra_%s%s_DDEC6.svg' % (cation[0],"all"), format="svg", dpi=2000, bbox_inches='tight')
plt.savefig('./figures_for_article/figure5_ChargeVSSpectra_%s%s_DDEC6.eps' % (cation[0],"all"), format="eps", dpi=2000, bbox_inches='tight')
