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

#Function for calculating correction factors for each atom
# filename - .xyz geometry file of ion pair
# DDEC - DDEC charges of ion pair
# RETURNS: Array of correction factors to be subtracted from binding energies of given ion pair
def ReadCoordinatesCalculateCorr(filename, DDEC):
    datafile = np.loadtxt(filename, skiprows=2, usecols=(1,2,3))
    x_coords = datafile[:,0]
    y_coords = datafile[:,1]
    z_coords = datafile[:,2]
    corrections = []
    for i in range(len (DDEC)):
        correction = 0
        for j in range(len (DDEC)):
            if(i!=j):
                correction+=float(DDEC[j])/(((x_coords[i]-x_coords[j])**2 + (y_coords[i]-y_coords[j])**2 + (z_coords[i]-z_coords[j])**2)**0.5)
#                print(DDEC[j])
#        print(correction)
        corrections+=[correction]
    return corrections

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
	

#Fitting
from scipy.optimize import curve_fit
def func(X, a, b,c): 
    chg, corr = X
    return a*chg+b*corr+c

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

DDEC_atom_long=[] #One long array of DDEC 6 charges to fit
DDEC_atom_long_array = [] #Array of DDEC 6 charges where each ion is in separate array

Corrections_long=[] #One long array of corrections to fit function
Corrections_long_array=[] #Array of corrections where each ion is in separate array

BEs_long=[] #One long array of BEs to fit function
BEs_long_array=[] #Array of binding energies where each ion is in separate array


count = -1
for cat in cation:
	for an in anion:
	    count += 1
	    atom, DDEC, bader, catan, DDEC_atom, bader_atom = [], [], [], [], [], []
	    atom, DDEC, bader, catan = ReadChargeFile('./data/' + cat + an + '/' + cat + an + '_charges.out')
	
	    BE_corrections = ReadCoordinatesCalculateCorr('./data/' + cat + an + '/' + cat + an + '.xyz', DDEC)
        #Cut list of corrections to leave only chosen atoms
	    BE_corrections_atom, idk = SortCharges(atom, BE_corrections, bader, whatAtom)

	
	    DDEC_cat, DDEC_an, bader_cat, bader_an = Charges(DDEC, bader, catan)
	    DDEC_atom, bader_atom = SortCharges(atom, DDEC, bader, whatAtom)
	    atomname, m = [], []
	    m, atomname = ReadSpectraFile('./data/' + cat + an + '/' + cat + an + '.out')
	    

        ####BE - correction
	    Corrections_long+=BE_corrections_atom
	    Corrections_long_array+=[BE_corrections_atom]
       
	    DDEC_atom_long+=DDEC_atom
	    DDEC_atom_long_array += [DDEC_atom]
        
	    BEs_long +=m[:len(DDEC_atom)]
	    BEs_long_array+=[m[:len(DDEC_atom)]]
		
	    for i in range(len(DDEC_atom)):
	    	nr_of_atom = i
	    	Y.append([DDEC_atom[nr_of_atom]])
	    	X.append([m[nr_of_atom]])
	    plt.scatter(DDEC_atom, m[:len(DDEC_atom)], c=ncolors[count], marker = marker[count], s=16, label=an)


##Fitting correction function and plotting calculated corrected values
X1=[np.asarray(DDEC_atom_long), np.asarray(Corrections_long)]
popt, pcov = curve_fit(func, X1, BEs_long)
#plt.plot(DDEC_atom_long,BEs_long, 'ko')
BEs_corr = func(X1, *popt)
#########Plot scattered points
count = -1
for cat in cation:
    for an in anion:
        count += 1
        X2=[np.asarray(DDEC_atom_long_array[count]), np.asarray(Corrections_long_array[count])]
        BEs_corr_short = func(X2, *popt)
        plt.scatter(np.asarray(DDEC_atom_long_array[count]), BEs_corr_short, c=ncolors[count], marker = "x", s=16, label=an)

from sklearn.linear_model import LinearRegression
###Linear correlation without correction
X = np.asarray(DDEC_atom_long)[:, np.newaxis]
reg1 = LinearRegression().fit(X, BEs_long) 
print("R2 "+ str(reg1.score(X, BEs_long)))
y = reg1.predict(X)
plt.plot(X, y, 'k', lw=1.0)
###Linear correlation with correction
reg2 = LinearRegression().fit(X, BEs_corr) 
print("R2 "+ str(reg2.score(X, BEs_corr)))
y2 = reg2.predict(X)
plt.plot(X, y2, 'k', lw=1.0)

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
