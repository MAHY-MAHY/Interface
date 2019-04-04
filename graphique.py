from math import*
from scipy import*
import numpy as np
import pickle
import csv
import matplotlib.pyplot as plt




fichier = np.loadtxt("outX", delimiter=' ')
fichier1 = np.loadtxt("outY", delimiter=' ')







plt.plot(fichier,fichier1)
plt.show()

