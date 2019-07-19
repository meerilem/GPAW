import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

numberlist = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
atomname, m = [], []
ncolors=['#fa5078', '#fa5078', '#fa5078', '#fa5078', '#fa5078', '#fa5078', '#fa5078', '#fa5078']

def ReadTheFile1(filename):
# notatom=True
    m, atomname = [], []
    with open(filename) as f:
        count = 0
        boolean = False
        for line in f:
            count +=1
            #This helps to jump over empty or random? values of atoms C-1 and N-1
            if boolean:
                boolean = False
                continue
#            if len(line.strip())==0:
#                continue
            tmp = line.split()
            if(tmp[0]=='C-1' or tmp[0]=='N-1'):
                boolean=True
                continue

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

#In dataframe are written the order number of carbon, which corresponds to alcylic carbon in cation
cationAlcylicCarbon = pd.DataFrame([[2,8,3,0,0]], columns = ['TEPA','BPy','Pyr14','BMIm','EMIm'])

anion = ['TFSI', 'FSI', 'PF6', 'BCN4', 'Cl', 'Br', 'I', 'BF4']
# marker = ["o", "v", "^", "p", "8", "*", "s", "d", "+"] 
marker = ["x", "x", "x", "x", "x", "x", "x", "x"] 
alpha = [0.45,0.45,0.45,0.45,0.45,0.45,0.45,0.45] 
#anion = [anion[1]]
#cat = cation[4]
namelist = []
#an = 'Cl'	
atom='C'
s=0.25
count = -1

ax = plt.figure(figsize=(8.3/2.54,8.3/2.54))

X, Y = [], []
for cat in cation:
    for an in anion:

        atomname, m = [], []
#        print(cat+an)
        m, atomname = ReadTheFile1('./data/' + cat + an + '/' + cat + an + '.out')	
        m_tmp, atomname = GetListOfAtoms(atom, m, atomname)
        m = [x for x in m_tmp] #Isnt' this redundant?? 
        
        #Shift the KS values such that aliphatic carbon is 285 eV
        shift = m[cationAlcylicCarbon[cat][0]] - 285 #eV
        m =(np.asarray(m) - shift).tolist()
        
        atomname1s, m1s = [], []
        m1s, atomname1s = ReadTheFile1('./data/' + cat + an + '/' + cat + an + '_1s.out')	
        m_tmp1s, atomname1s = GetListOfAtoms(atom, m1s, atomname1s)
        m1s = [x1 for x1 in m_tmp1s]
        
        #Shift the 1s values such that aliphatic carbon is 285 eV
        shift = m1s[cationAlcylicCarbon[cat][0]] - 285 #eV
        m1s =(np.asarray(m1s) - shift).tolist()
        
        for i in range(len(m1s)):
            nr_of_atom = i
            Y.append([m[nr_of_atom]]) #Here are deltaKS
            X.append([m1s[nr_of_atom]]) #Here are 1s
        if an == 'BF4' and cat == 'EMIm':
            bf4 = plt.scatter(m, m1s, c='#7fffbb', marker = "D", s=24, label=an,alpha=1,edgecolors="#333333",linewidths=0.7)
        else:
            plt.scatter(m, m1s, c=ncolors[count], marker = marker[count], s=16, label=an,alpha=alpha[count],edgecolors="#333333",linewidths=0.7)
    count += 1

    
from sklearn.linear_model import LinearRegression
reg = LinearRegression().fit(X, Y) 
xplot = np.linspace(start=282,stop=294,num=100).reshape(-1,1)
y = reg.predict(xplot)
plt.plot(y, xplot, 'k-', lw=1.0,linestyle='--')
params = {'legend.fontsize': 7}
plt.rcParams.update(params)
print("R2 value = " + str(reg.score(X,Y)))

y_est = reg.predict(X)

differenceSquared = np.square(np.asarray(y_est)-np.asarray(Y))
std_error = np.sqrt(sum(differenceSquared)/len(differenceSquared))
print('The standard error of the estimate '+ str(std_error[0]))

# from matplotlib.legend import Legend

# leg = Legend(ax,[bf4],[r"BF$_{4}$"],loc=(0.22,0.85),
# 	frameon=False,fancybox=True, fontsize=7,handlelength=0)
# ax.add_artist(leg)

plt.xlim(left=282, right=290)
plt.ylim(bottom=282, top=290)

plt.xticks(fontsize=7)
plt.yticks(fontsize=7)
plt.yticks([282,284,286,288,290])
plt.minorticks_on()
plt.xlabel(r"$\Delta$" + "KS binding energy / eV", fontsize=8)
plt.ylabel("1s binding energy / eV", fontsize=8)
plt.tight_layout()
plt.savefig('./figures_for_article/figure4.png', dpi=600)
plt.savefig('./figures_for_article/figure4.svg', dpi=600)
#plt.savefig('./figures_for_article/figure6_Spectra1sVSSpectradKs.eps'), dpi=2000, bbox_inches='tight')