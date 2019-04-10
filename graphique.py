import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
from pylab import*
from matplotlib.path import Path
from matplotlib.patches import PathPatch

#coor=np.loadtxt("outX.txt");
TMP=np.loadtxt("outY.txt");
TMP2=np.loadtxt("outYex.txt");
print("entrez le nombre de mailles en X")
i=input();
print("entrez le temps final")
t=input()
taillef=int(i/(t+1))
#while (coor[0,1]==coor[i,1]):
#    i=i+1;
A=np.zeros((taillef,taillef));
C=np.zeros((taillef,taillef))
B=np.zeros((taillef,taillef));
for j in range(taillef):
    for k in range(taillef):
        A[j,k]=TMP[i*i-1-k-i*j]*(1-TMP[i*i-1-k-j*i])+0
        C[j,k]=TMP[i*i-1-k-j*i]+0
        B[j,k]=TMP2[i*i-1-k-j*i]+0

x=np.loadtxt("outX.txt")
X,Y = np.meshgrid(x[0:taillef,0],x[0:taillef,0])
plt.pcolor(X, Y, C)
plt.colorbar()
plt.show()

plt.pcolor(X, Y, B)
plt.colorbar()
plt.show()

plt.pcolor(X, Y, A)
plt.colorbar()
plt.show()
