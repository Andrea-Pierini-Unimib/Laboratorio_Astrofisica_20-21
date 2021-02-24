import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from timeit import default_timer as timer
import datetime
from astropy.cosmology import LambdaCDM
import astropy.units as u
from sys import stdout

def lum_dist(z):
    obj = LambdaCDM(H0=69.6, Om0=0.286, Ode0=1-0.286)
    dist = obj.luminosity_distance(z).value * 3.08567758e24 # cm
    return dist

def lum_flux(f,z):
    return f * 4 * np.pi * lum_dist(z)**2 # erg s^-1 Hz^-1

def lum_dist_MPC(z):
    obj = LambdaCDM(H0=69.6, Om0=0.286, Ode0=1-0.286)
    dist = obj.luminosity_distance(z).value # MPC
    return dist

def com_dist_MPC(z):
    obj = LambdaCDM(H0=69.6, Om0=0.286, Ode0=1-0.286)
    dist = obj.comoving_distance(z).value # MPC
    return dist

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

#####################################################################################

N = len(ID)
f_FUV = np.zeros(N)
f_NUV = np.zeros(N)
FUV_sx = 1344
FUV_dx = 1786
NUV_sx = 1771
NUV_dx = 2831
j = 0
for a_file in ID:
    fitsfile = fits.open(a_file)
    data = fitsfile[1].data
    L = data.shape[0]
    wave_rest = np.zeros(L)
    f_lambda_rest = np.zeros(L)
    for i in range(L):
        wave_rest[i] = data[i][1]
        f_lambda_rest[i] = data[i][3]
    wave_rest = np.array(wave_rest)
    f_lambda_rest = np.array(f_lambda_rest)
    f_fuv = f_lambda_rest[(wave_rest >= 1344) & (wave_rest <= 1786)]
    f_FUV[j] = np.mean(f_fuv)
    f_nuv = f_lambda_rest[(wave_rest >= 1771) & (wave_rest <= 2831)]
    f_NUV[j] = np.mean(f_nuv)
    j += 1
    stdout.write('\r Progress = {:.3f}%'.format(j/N*100.)) # Dynamic print
    fitsfile.close()

#####################################################################################

# cerchiamo comoving volume di astropy

#####################################################################################

# stiamo assumendo che f_NUV sia in: erg cm^-2 s^-1 Å^-1

c = 2.99792458e18 # Å / s
lambda_nuv = 2271
lambda_fuv = 1528

fl_FUV = f_FUV * lambda_fuv**2 / c
fl_NUV = f_NUV * lambda_nuv**2 / c

# ora fl_FUV è in: erg / cm^-2 s^-1 Hz^-1

L_FUV = np.zeros(len(Z))
L_NUV = np.zeros(len(Z))

for i in range(len(Z)):
    L_FUV[i] = lum_flux(fl_FUV[i], Z[i])
    L_NUV[i] = lum_flux(fl_NUV[i], Z[i])

# dai flussi [cgs] siamo passati alle luminosità [cgs]

fig , (ax_fuv , ax_nuv) = plt.subplots(2,1,figsize=(15,15))

logbins = np.logspace(25,30,11)

ax_fuv.hist(L_FUV, bins=logbins, edgecolor='black', color='lightgrey')
ax_nuv.hist(L_NUV, bins=logbins, edgecolor='black', color='lightgrey')

ax_fuv.set_xscale('log')
ax_nuv.set_xscale('log')
ax_fuv.set_yscale('log')
ax_nuv.set_yscale('log')
ax_fuv.grid(ls=':',which='both')  
ax_nuv.grid(ls=':',which='both')
ax_fuv.set_xlabel('$L_{FUV}$ [erg s$^{-1}$ Hz$^{-1}$]', fontsize=15)
ax_nuv.set_xlabel('$L_{NUV}$ [erg s$^{-1}$ Hz$^{-1}$]', fontsize=15)
ax_fuv.set_ylabel('$counts$', fontsize=15)
ax_nuv.set_ylabel('$counts$', fontsize=15)

fig.savefig('PNG/0_ALL_FIELDS_histo_luminosity.png', dpi=600)

#####################################################################################

# fissata la galassia: Z[i] , L_FUV[i] , L_NUV[i] , M[:,i] , SFR[:,i] , ID_f[i]
# ordinare per redshift:
    
sort_index = np.argsort(Z)

Z     = Z[sort_index]
L_FUV = L_FUV[sort_index]
L_NUV = L_NUV[sort_index]
# ID_f  = ID_f[sort_index]
ID    = ID[sort_index]

for i in range(3):
    M[i]   = M[i][sort_index]
    SFR[i] = SFR[i][sort_index]

str_Z     = []
str_L_FUV = []
str_L_NUV = []
str_ID_f = []
str_ID = []
str_M = [[] for i in range(3)]
str_SFR = [[] for i in range(3)]
bin_1 = (Z <= 1)
bin_2 = (Z > 1) & (Z <= 2)
bin_3 = (Z > 2)
bins  = [bin_1, bin_2, bin_3]

for i in range(len(bins)):
    str_Z.append(Z[bins[i]])
    str_L_FUV.append(L_FUV[bins[i]])
    str_L_NUV.append(L_NUV[bins[i]])
    str_ID_f.append(ID_f[bins[i]])
    str_ID.append(ID[bins[i]])
    for j in range(3):
        str_M[j].append(M[j][bins[i]])
        str_SFR[j].append(SFR[j][bins[i]])

#####################################################################################

# controllare il lim_inf del primo bin

d_omega = np.array([3.392190201880897e-05, 
           5.251153102188274e-07, 
           3.559951792354943e-05, 
           1.1785411038838987e-05, 
           1.5116216487179068e-06
          ])
D_omega = np.sum(d_omega)
# controllare nei paper l'area in °^2 

NZ = 3
d_RA_bins = np.zeros(NZ)
d_DEC_bins = np.zeros(NZ)
Z_bins = np.array([0,1,2,4])
d_Z_bins = lum_dist_MPC(Z_bins)
d_Z_bins_com = com_dist_MPC(Z_bins)

# volumi dei singoli tronchi di piramide 
V_1 = 1./3. * d_Z_bins_com[1]**3 * D_omega
V_2 = 1./3. * d_Z_bins_com[2]**3 * D_omega - V_1
V_3 = 1./3. * d_Z_bins_com[3]**3 * D_omega - V_1 - V_2

# 1./3. * d_z * ra * d_z * dec * sin(dec) * d_z
# dV = r^2 dr domega
# V = integro in dr = 1/3 r^3 domega

print('V_1:', V_1)
print('V_2:', V_2)
print('V_3:', V_3)

V = np.array([V_1, V_2, V_3])

print('Total number of galaxies = {:d}'.format(len(ID)))

######################################################
######################################################
## LA COSA PIU' SERIA: LUMINOSITY FUNCTION (Z BINS) ##
######################################################
######################################################

# mettere errori possioniani

L_sun = 3.826e33 # erg s^-1
# L_sun_FUV # erg s^-1 Hz^-1
# L_sun_NUV # erg s^-1 Hz^-1
# sostituire con la luminosità del sole nella stessa banda 
# oppure passare alle magnitudini assolute

fig , ax = plt.subplots(1, 2, figsize=(12,6), constrained_layout=True)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"]
})

fig.suptitle('Luminosity function - redshift bins', fontsize=15)
ax[0].set_title('Far-UV band', fontsize=15)
ax[1].set_title('Near-UV band', fontsize=15)

logbins = np.logspace(25,30,16) / L_sun   # np.logspace(25,30,11)
centres = np.logspace(25,30,31) / L_sun   # np.logspace(25,30,21)
centres = centres[1:30:2]                 # centres[1:20:2]
widths  = np.zeros(15)
for i in range(len(logbins)-1):
    widths[i] = logbins[i+1] - logbins[i]

density_FUV = []
density_NUV  = []
colorstring = ['red','orange','yellow']
labelstring = ['$0 < z \leq 1$','$1 < z \leq 2$','$2 < z \leq 4$']
for i in range(3):
    histo_fuv, trash_1 = np.histogram(str_L_FUV[i] / L_sun, bins=logbins)
    histo_nuv, trash_2 = np.histogram(str_L_NUV[i] / L_sun, bins=logbins)
    for j in range(len(histo_fuv)):
        histo_fuv[j] = histo_fuv[j] # / widths[j]
        histo_nuv[j] = histo_nuv[j] # / widths[j]
    density_FUV.append(histo_fuv / V[i])
    density_NUV.append(histo_nuv / V[i])
    fuv_index = density_FUV[i] != 0.
    nuv_index = density_NUV[i] != 0.
    fuv_centres = centres[fuv_index]
    nuv_centres = centres[nuv_index]
    fuv_test = density_FUV[i][fuv_index]
    nuv_test = density_NUV[i][nuv_index]
    ax[0].plot(fuv_centres, fuv_test, color=colorstring[i])
    ax[0].scatter(centres, density_FUV[i], color=colorstring[i], label=labelstring[i])
    ax[1].plot(nuv_centres, nuv_test, color=colorstring[i])
    ax[1].scatter(centres, density_NUV[i], color=colorstring[i], label=labelstring[i])


for j in range(2):
    ax[j].set_xscale('log')
    ax[j].set_yscale('log')
    ax[j].set_xlim(1e-9,1e-3)
    # ax[j].set_ylim(1e-4,1e3)
    ax[j].set_ylim(1e-7,1e-2)
    ax[j].grid(ls=':',which='both')
    ax[j].set_ylabel('$n$ [Mpc$^{-3}$ ($L_{\odot}$ Hz$^{-1}$)$^{-1}$]', fontsize=15) 
    ax[j].legend(frameon=True, loc=2)
ax[0].set_xlabel('$L_{FUV}$ [$L_{\odot}$ Hz$^{-1}$]', fontsize=15)
ax[1].set_xlabel('$L_{NUV}$ [$L_{\odot}$ Hz$^{-1}$]', fontsize=15)

fig.savefig('PNG/0_ALL_FIELDS_luminosity_function_z_bins.png', dpi=600)
fig.savefig('EPS/0_ALL_FIELDS_luminosity_function_z_bins.eps')
fig.savefig('PDF/0_ALL_FIELDS_luminosity_function_z_bins.pdf')

#####################################################################################

ssfr = np.zeros_like(SFR)

for j in range(3):
    for i in range(len(SFR[1])):
        ssfr[j][i] = 10**SFR[j][i] / (10**M[j][i])

ssfr_bins = np.logspace(-10,-8,3)
ssfr_lim_sx = np.amin(ssfr[1])
ssfr_trubins = np.zeros(len(ssfr_bins)+1)
ssfr_trubins[0] = ssfr_lim_sx
for i in range(len(ssfr_bins)):
    ssfr_trubins[1+i] = ssfr_bins[i]
    
sfr_index = np.argsort(ssfr[1])

Z     = Z[sfr_index]
L_FUV = L_FUV[sfr_index]
L_NUV = L_NUV[sfr_index]
ID_f  = ID_f[sfr_index]
ID    = ID[sfr_index]

for i in range(3):
    ssfr[i] = ssfr[i][sfr_index]

sfr_Z     = []
sfr_L_FUV = []
sfr_L_NUV = []
sfr_ID_f = []
sfr_ID = []
sfr_M = []
sfr_SSFR = [[] for i in range(3)]

sbin_1 = (ssfr[1] <= ssfr_trubins[1])
sbin_2 = (ssfr[1] > ssfr_trubins[1]) & (ssfr[1] <= ssfr_trubins[2])
sbin_3 = (ssfr[1] > ssfr_trubins[2])
sbins  = [sbin_1, sbin_2, sbin_3]

for i in range(len(sbins)):
    sfr_Z.append(Z[sbins[i]])
    sfr_L_FUV.append(L_FUV[sbins[i]])
    sfr_L_NUV.append(L_NUV[sbins[i]])
    sfr_ID_f.append(ID_f[sbins[i]])
    sfr_ID.append(ID[sbins[i]])
    for j in range(3):
        sfr_SSFR[j].append(ssfr[j][sbins[i]])

#########################################################
#########################################################
## LA COSA PIU' SERIA: LUMINOSITY FUNCTION (SSFR BINS) ##
#########################################################
#########################################################


L_sun = 3.826e33 # erg/s


fig , ax = plt.subplots(1, 2, figsize=(12,6), constrained_layout=True)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"]
})

fig.suptitle('Luminosity function - specific star formation rate bins', fontsize=15)
ax[0].set_title('Far-UV band', fontsize=15)
ax[1].set_title('Near-UV band', fontsize=15)

logbins = np.logspace(25,30,16) / L_sun     # np.logspace(25,30,11)
centres = np.logspace(25,30,31) / L_sun     # np.logspace(25,30,21)
centres = centres[1:30:2]                   # centres[1:30:2]
widths  = np.zeros(15)                      # np.zeros(10)
for i in range(len(logbins)-1):
    widths[i] = logbins[i+1] - logbins[i]

V_tot = V[0] + V[1] + V[2]

sfr_density_FUV = []
sfr_density_NUV  = []
colorstring = ['midnightblue','dodgerblue','lightskyblue']
labelstring = ['$10^{-16} < ssfr \leq 10^{-10}$','$10^{-10} < ssfr \leq 10^{-9}$','$10^{-9} < ssfr \leq 10^{-8}$']
for i in range(3):
    histo_fuv, trash_1 = np.histogram(sfr_L_FUV[i] / L_sun, bins=logbins)
    histo_nuv, trash_2 = np.histogram(sfr_L_NUV[i] / L_sun, bins=logbins)
    for j in range(len(histo_fuv)):
        histo_fuv[j] = histo_fuv[j] # / widths[j]
        histo_nuv[j] = histo_nuv[j] # / widths[j]
    sfr_density_FUV.append(histo_fuv / V_tot)
    sfr_density_NUV.append(histo_nuv / V_tot)
    fuv_index = sfr_density_FUV[i] != 0.
    nuv_index = sfr_density_NUV[i] != 0.
    fuv_centres = centres[fuv_index]
    nuv_centres = centres[nuv_index]
    fuv_test = sfr_density_FUV[i][fuv_index]
    nuv_test = sfr_density_NUV[i][nuv_index]
    ax[0].plot(fuv_centres, fuv_test, color=colorstring[i])
    ax[1].plot(nuv_centres, nuv_test, color=colorstring[i])
    ax[0].scatter(centres, sfr_density_FUV[i], color=colorstring[i], label=labelstring[i])
    ax[1].scatter(centres, sfr_density_NUV[i], color=colorstring[i], label=labelstring[i])

for j in range(2):
    ax[j].set_xscale('log')
    ax[j].set_yscale('log')
    ax[j].set_xlim(1e-9,1e-3)
    # ax[j].set_ylim(1e-4,1e3)
    ax[j].set_ylim(1e-7,2e-4)
    ax[j].grid(ls=':',which='both')
    ax[j].set_ylabel('$n$ [Mpc$^{-3}$ ($L_{\odot}$ Hz$^{-1}$)$^{-1}$]', fontsize=15) 
    ax[j].legend(frameon=True)
ax[0].set_xlabel('$L_{FUV}$ [$L_{\odot}$ Hz$^{-1}$]', fontsize=15)
ax[1].set_xlabel('$L_{NUV}$ [$L_{\odot}$ Hz$^{-1}$]', fontsize=15)

fig.savefig('PNG/0_ALL_FIELDS_luminosity_function_ssfr_bins.png', dpi=600)
fig.savefig('EPS/0_ALL_FIELDS_luminosity_function_ssfr_bins.eps')
fig.savefig('PDF/0_ALL_FIELDS_luminosity_function_ssfr_bins.pdf')