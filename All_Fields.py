import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from timeit import default_timer as timer
import datetime
from astropy.cosmology import LambdaCDM
import astropy.units as u
from sys import stdout

print("Hello Twats")

#####################################################################################

listname = ['/data1/astrolab/data/projects/3dhst_20_21/SPS/AEGIS/Plots/Fit_values.txt',
           '/data1/astrolab/data/projects/3dhst_20_21/SPS/COSMOS/Plots/Fit_values.txt',
           '/data1/astrolab/data/projects/3dhst_20_21/SPS/GOODSN/Plots/Fit_values.txt',
           '/data1/astrolab/data/projects/3dhst_20_21/SPS/GOODSS/Plots/Fit_values.txt',
           '/data1/astrolab/data/projects/3dhst_20_21/SPS/UDS/Plots/Fit_values.txt']
datapath = ['/data1/astrolab/data/projects/3dhst_20_21/SPS/AEGIS/Plots/id',
           '/data1/astrolab/data/projects/3dhst_20_21/SPS/COSMOS/Plots/id',
           '/data1/astrolab/data/projects/3dhst_20_21/SPS/GOODSN/Plots/id',
           '/data1/astrolab/data/projects/3dhst_20_21/SPS/GOODSS/Plots/id',
           '/data1/astrolab/data/projects/3dhst_20_21/SPS/UDS/Plots/id']
extension = '_best_sed.fits'
ID = []
ID_f = []
Z = []
SFR = [[] for i in range(3)]
M = [[] for i in range(3)]

for i in range(len(listname)):
    name = listname[i]
    path = datapath[i]
    listfile = open(name, 'r')
    line_0 = listfile.readline()
    for line in listfile:
        a = line.strip('\n')
        a = a.split()
        ID.append(path+a[0]+extension)
        ID_f.append(float(a[0]))
        Z.append(float(a[1]))
        for j in range(3):
            SFR[j].append(float(a[23+j]))
            M[j].append(float(a[14+j]))
    listfile.close()

ID = np.array(ID)
ID_f = np.array(ID_f)
Z = np.array(Z)
SFR = np.array(SFR)
M = np.array(M)