
# coding: utf-8

# # Figure 5 Binding energies

# In[1]:


from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


numberlist = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
atomname, m = [], []
ncolors=['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628', '#984ea3', '#999999', '#e41a1c', '#dede00']
#ncolors=['#CC3333', '#CC3333', '#CC3333', '#CC3333', '#CC3333', '#CC3333', '#CC3333', '#CC3333']

polarizabilities={'H':4.50711,'C': 11.3, 'N': 7.4, 'O':5.3, 'F':3.74, 'B':20.5, 'S':19.4,'P':25,'Cl':14.6, 'Br':21, 'I':32.9}

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

#Function for calculating correction factors for each atom from their neigbours charges
# coordinates - 2D Array of coordinates of ion pair [atom_name,x,y,z]
# DDEC - DDEC charges of ion pair
# dipoles - 2D Array of x,y,z atomic dipoles
# RETURNS: Array of correction factors calculated from neighbors charge
def CalculateCorrCharge(coordinates, DDEC):

    atom_name=coordinates[0]
    x_coords = coordinates[1]
    y_coords = coordinates[2]
    z_coords = coordinates[3]
    
    corrections = []
    for i in range(len (DDEC)):
        correction = 0
        for j in range(len (DDEC)):
            if(i!=j):
#                print(polariza[atom_name[i]])
#                correction+=float(DDEC[j])*polariza[atom_name[j]]/(((x_coords[i]-x_coords[j])**2 + (y_coords[i]-y_coords[j])**2 + (z_coords[i]-z_coords[j])**2)**0.5)
                correction+=float(DDEC[j])/(((x_coords[i]-x_coords[j])**2 + (y_coords[i]-y_coords[j])**2 + (z_coords[i]-z_coords[j])**2)**0.5)

#                print(DDEC[j])
        

        corrections+=[correction]
#    print(correction)
    return corrections

#Function for reading in coordinates
# filename - .xyz geometry file of ion pair
# RETURNS: Array of [atom_name, x-coordinates, y-coordinates, z-coordinates]
def ReadCoordinates(filename):
    datafile = np.loadtxt(filename, skiprows=2, usecols=(1,2,3))
    atom_name=np.loadtxt(filename, skiprows=2,usecols=(0), dtype=str)
    x_coords = datafile[:,0]
    y_coords = datafile[:,1]
    z_coords = datafile[:,2]
    
    return [atom_name, x_coords, y_coords, z_coords]

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


# In[2]:


#Fitting only charge
def func1(Z, c): 
    chg, corr = Z
    return 18.8*chg+c


# In[3]:


#Fitting only charge and neighbours charge
def func2(Z, c): 
    chg, corr = Z
    return 18.8*chg+14.4*corr+c


# In[4]:


#Fitting charge, neighbours charge and dipoles
from scipy.optimize import curve_fit
def func3(Z,c): 
    chg, corr1, corr2 = Z
    return 18.8*chg+14.4*corr1+14.4*corr2+c


# In[5]:


#Function to read in dipoles
#filename - DDEC6 charge file (shortened)
#returns dipoles in 2D array, units are e*A

def readDipoles(filename):
    
    au2Cm = 8.4783536255*10**-30
    elChg=1.60217662*10**-19 #C/e
    m2A = 10**10 #A/m
    au2eA=au2Cm*m2A/elChg
    
    columns = np.loadtxt(filename, skiprows=2, usecols=(6,7,8))
    return [list(columns[:,0]*au2eA), list(columns[:,1]*au2eA), list(columns[:,2]*au2eA)]


# In[6]:


#Function for calculating correction factors for each atom from their neigbours charges
# coordinates - 2D Array of coordinates of ion pair [atom_name,x,y,z]
# DDEC - DDEC charges of ion pair
# dipoles - 2D Array of x,y,z atomic dipoles
# RETURNS: Array of correction factors calculated from neighbors charge
def CalculateCorrDip(coordinates, DDEC, dipoles):

    atom_name=coordinates[0]
    x_coords = coordinates[1]
    y_coords = coordinates[2]
    z_coords = coordinates[3]
    
    dip_x = dipoles[0]
    dip_y = dipoles[1]
    dip_z = dipoles[2]
    
    corrections = []
    for i in range(len (DDEC)):
        correction = 0
        for j in range(len (DDEC)):
            if(i!=j):
                correction+=float(DDEC[i])*((x_coords[i]-x_coords[j])*dip_x[j]+(y_coords[i]-y_coords[j])*dip_y[j]+(z_coords[i]-z_coords[j])*dip_z[j])/(((x_coords[i]-x_coords[j])**2 + (y_coords[i]-y_coords[j])**2 + (z_coords[i]-z_coords[j])**2)**(3*0.5))
        corrections+=[correction]

    return corrections


# # BE vs f(q, d)

# In[7]:


f1, (ax) = plt.subplots(1, 1, sharey=False,sharex=True, figsize=(8.3/2.54, 8.3/2.54))
f1.subplots_adjust(hspace=0)


# In[8]:


df=pd.DataFrame(columns=['Cation', 'Anion', 'q R2', 'f(q) R2', 'f(q,d) R2'])


# In[9]:


ncolors=['#aa0a3c', '#aa0a3c', '#aa0a3c', '#aa0a3c', '#aa0a3c', '#aa0a3c', '#aa0a3c', '#aa0a3c']
#ncolors=["#a0fa82"]*8

ncolors=["#33cc33"]*8

########################################################
cation = ['TEPA', 'BPy', 'Pyr14', 'BMIm', 'EMIm'] #Cations array
anion = ['TFSI', 'FSI', 'BF4', 'PF6', 'BCN4', 'Cl', 'Br', 'I'] #Anions Array

marker = ["x", "x", "x", "x", "x", "x", "x", "x"] 
alpha = [0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45] 

#anion = ['TFSI', 'FSI', 'BF4', 'PF6', 'BCN4']
#anion = [anion[2]]
cation = [cation[4]]

X, Y = [], []
nr_of_atom = 4
whatAtom = 'C'


DDEC_atom_long=[] #One long array of DDEC 6 charges to fit
DDEC_atom_long_array = [] #Array of DDEC 6 charges where each ion is in separate array

Corrections_longCharge=[] #One long array of corrections to fit function
Corrections_longCharge_array=[] #Array of corrections where each ion is in separate array

Corrections_longDipole=[] #One long array of corrections to fit function
Corrections_longDipole_array=[] #Array of corrections where each ion is in separate array

BEs_long=[] #One long array of BEs to fit function
BEs_long_array=[] #Array of binding energies where each ion is in separate array


specialAnion="BF4" #The anion we want to show

count = -1
#Iterating over list of cations and anions
for cat in cation:

    for an in anion:

        count += 1
        
        #Create empty arrays and load in data from files 
        ##Atom - atom names
        ##DDEC - atomic charges by DDEC
        ##bader - atomic charges by Bader
        ##catan - AN/CAT
        
        atom, DDEC, bader, catan, DDEC_atom, bader_atom = [], [], [], [], [], []

        atom, DDEC, bader, catan = ReadChargeFile('./data/' + cat + an + '/' + cat + an + '_charges.out')
        
        #Read in dipoles from DDEC .xyz (only returns dipoles)
        dipoles = readDipoles('./data/' + cat + an + '/DDEC6_charge.xyz')
        
        #Read in coordinates from .xyz file and calculate correction
        coords = ReadCoordinates('./data/' + cat + an + '/' + cat + an + '.xyz')
        
        #Calculate corrections from coordinates, charge and dipoles
        BE_correctionsCharge = CalculateCorrCharge(coords, DDEC) #Correction of neighboring atoms charges
        BE_correctionsDipole = CalculateCorrDip(coords, DDEC, dipoles) #Correction of neighboring atoms dipoles
        
        #Cut list of corrections to leave only chosen atoms
        BE_corrections_atomCharge, idk = SortCharges(atom, BE_correctionsCharge, bader, whatAtom)
        BE_corrections_atomDipole, idk = SortCharges(atom, BE_correctionsDipole, bader, whatAtom)

        
        DDEC_cat, DDEC_an, bader_cat, bader_an = Charges(DDEC, bader, catan)
        DDEC_atom, bader_atom = SortCharges(atom, DDEC, bader, whatAtom)

        atomname, m = [], []

        m, atomname = ReadSpectraFile('./data/' + cat + an + '/' + cat + an + '.out')

        ####BE - correction

        Corrections_longCharge+=BE_corrections_atomCharge
        Corrections_longCharge_array+=[BE_corrections_atomCharge]

        Corrections_longDipole+=BE_corrections_atomDipole
        Corrections_longDipole_array+=[BE_corrections_atomDipole]

        
        DDEC_atom_long+=DDEC_atom
        DDEC_atom_long_array += [DDEC_atom]
        
        BEs_long +=m[:len(DDEC_atom)]
        BEs_long_array+=[m[:len(DDEC_atom)]]


# In[10]:


#########Plot scattered points
count = -1
for cat in cation:
    for an in anion:
        count += 1
        X2=[np.asarray(DDEC_atom_long_array[count]), np.asarray(Corrections_longCharge_array[count])]
        ###Linear correlation
        X = np.asarray(BEs_long_array[count])[:,np.newaxis]#np.asarray(BEs_long_array[count])[:, np.newaxis]
        x_range=np.arange(289,296,0.1)

        print(an)
        
        ##Correlation of ion pair with different corrections
        
        ##Only charge
        popt_q, pcov_q = curve_fit(func1, X2, BEs_long_array[count])
        BEs_corr_charge = func1(X2, *popt_q)
        reg_q = LinearRegression().fit(X, BEs_corr_charge)
        y_q=reg_q.predict(x_range[:, np.newaxis])
        print("q R2=" + str(reg_q.score(X, BEs_corr_charge)))       

        
        ##Charge of surrounding atoms
        popt_qq, pcov_qq = curve_fit(func2, X2, BEs_long_array[count])
        BEs_corr_charge_q = func2(X2, *popt_qq)
        reg_qq = LinearRegression().fit(X, BEs_corr_charge_q)
        y_qq=reg_qq.predict(x_range[:, np.newaxis])        
        print("f(q) R2=" + str(reg_qq.score(X, BEs_corr_charge_q)))      

        ##Charge, neighbours, dipoles
        X3 = [np.asarray(DDEC_atom_long_array[count]), np.asarray(Corrections_longCharge_array[count]),
              np.asarray(Corrections_longDipole_array[count])]
        popt_qqd, pcov_qqd = curve_fit(func3, X3, BEs_long_array[count])
        BEs_corr_charge_qd = func3(X3, *popt_qqd)
        reg_qqd = LinearRegression().fit(X, BEs_corr_charge_qd)
        y_qqd=reg_qqd.predict(x_range[:, np.newaxis])           
        print("f(q,d) R2=" + str(reg_qqd.score(X, BEs_corr_charge_qd)))
        
        #Saving correlation factors to dataframe
        df2 = pd.DataFrame([[cat, an, reg_q.score(X, BEs_corr_charge), reg_qq.score(X, BEs_corr_charge_q), reg_qqd.score(X, BEs_corr_charge_qd)]],
                           columns=['Cation', 'Anion', 'q R2', 'f(q) R2', 'f(q,d) R2'])
        df=df.append(df2, ignore_index=True)

        if(an==specialAnion):


            ax.scatter(BEs_long_array[count], BEs_corr_charge,
                       color='#aa0a3c',
                       marker='D',alpha=1, s=16,
                       edgecolors="#333366",linewidths=0.7,label="q")
#            ax.scatter(BEs_long_array[count], BEs_corr_charge_q,
#                       color='#4daf4a',
#                       marker='D',alpha=1, s=16,edgecolors="#333366",linewidths=0.7,
#                       label="q+q_i")
#            ax.plot(x_range, y_qq, color='#4daf4a')

            ax.scatter(BEs_long_array[count], BEs_corr_charge_qd,
                       color='#4daf4a',
                       marker='D',alpha=1, s=16,edgecolors="#333366",linewidths=0.7,
                       label="q+q_i+dip")
            ax.plot(x_range, y_qqd, color='#4daf4a')

            

            #Create a DataFrame for checking data
#            df['q'] = DDEC_atom_long_array[count]
#            df['a*q'] = 18.8*np.asarray(DDEC_atom_long_array[count])
#            df['sum(q_i/r)'] = Corrections_longCharge_array[count]
#            df['b*sum(q_i/r)'] = 14.4*np.asarray(Corrections_longCharge_array[count])
#            df['Estimated BE charge'] = BEs_corr_charge
#            df['Estimated BE charge + neighbours'] = BEs_corr_charge_q
#            df['Delta KS'] = BEs_long_array[count]
            
        else:
            #Not corrected points
            ax.scatter(BEs_long_array[count],#BEs_long_array[count],
                        BEs_corr_charge, c='#aa0a3c', marker = marker[count],
                    alpha=alpha[count], s=16, label='_nolegend_',
                       edgecolors="#333366",
                       linewidths=0.7)
            
            #Corrected points
            ax.scatter(BEs_long_array[count],#BEs_long_array[count],
                        BEs_corr_charge_q, c='#4daf4a',#"#33cc33",
                       marker = marker[count],
                    alpha=alpha[count], s=16, label='_nolegend_',
                       edgecolors="#333366",
                       linewidths=0.7)



ax.tick_params(axis='x', which='both', direction='inout')
ax.set_xlabel(r"$\Delta$" + "KS binding energy / eV",fontsize=8)
ax.set_ylabel(r"$f(q,d)$ BE / eV",fontsize=8)
#ax.legend(prop={'size': 8})


# In[11]:


f1


# In[12]:


##Correlation of Set with different corrections

print("Set")

##Only charge
popt_q, pcov_q = curve_fit(func1, [np.asarray(DDEC_atom_long), np.asarray(Corrections_longCharge)],
                           BEs_long)
BEs_corr_charge = func1([np.asarray(DDEC_atom_long), np.asarray(Corrections_longCharge)], *popt_q)
reg_q = LinearRegression().fit(np.asarray(BEs_long)[:,np.newaxis], BEs_corr_charge)
y_q=reg_q.predict(x_range[:, np.newaxis])
print("q R2=" + str(reg_q.score(np.asarray(BEs_long)[:,np.newaxis], BEs_corr_charge)))      

        
##Charge of surrounding atoms
popt_qq, pcov_qq = curve_fit(func2, [np.asarray(DDEC_atom_long), np.asarray(Corrections_longCharge)], BEs_long)
BEs_corr_charge_q = func2([np.asarray(DDEC_atom_long), np.asarray(Corrections_longCharge)], *popt_qq)
reg_qq = LinearRegression().fit(np.asarray(BEs_long)[:,np.newaxis], BEs_corr_charge_q)
y_qq=reg_qq.predict(x_range[:, np.newaxis])        
print("f(q) R2=" + str(reg_qq.score(np.asarray(BEs_long)[:,np.newaxis], BEs_corr_charge_q)))       

##Charge, neighbours, dipoles
S3 = [np.asarray(DDEC_atom_long), np.asarray(Corrections_longCharge),
        np.asarray(Corrections_longDipole)]
popt_qqd, pcov_qqd = curve_fit(func3, S3, BEs_long)
BEs_corr_charge_qd = func3(S3, *popt_qqd)
reg_qqd = LinearRegression().fit(np.asarray(BEs_long)[:,np.newaxis], BEs_corr_charge_qd)
y_qqd=reg_qqd.predict(x_range[:, np.newaxis])           
print("f(q,d) R2=" + str(reg_qqd.score(np.asarray(BEs_long)[:,np.newaxis], BEs_corr_charge_qd)))


# In[14]:


df2 = pd.DataFrame([["Set", "Set", reg_q.score(np.asarray(BEs_long)[:,np.newaxis], BEs_corr_charge), reg_qq.score(np.asarray(BEs_long)[:,np.newaxis], BEs_corr_charge_q), reg_qqd.score(np.asarray(BEs_long)[:,np.newaxis], BEs_corr_charge_qd)]],
                    columns=['Cation', 'Anion', 'q R2', 'f(q) R2', 'f(q,d) R2'])
df=df.append(df2, ignore_index=True)
df.to_excel("./CorrelationData.xlsx")


# In[15]:


#Save figure
f1.savefig('./figures_for_article/figure5_ChargeVSSpectra_%s%s_DDEC6.png' % (cation[0],"all"), format="png", dpi=300, bbox_inches='tight')
f1.savefig('./figures_for_article/figure5_ChargeVSSpectra_%s%s_DDEC6.png' % (cation[0],"all"), format="png", dpi=300, bbox_inches='tight')

