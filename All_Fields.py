import numpy             as np
import matplotlib.pyplot as plt
from   astropy.io        import fits
from   timeit            import default_timer as timer
import datetime
from   astropy.cosmology import LambdaCDM
import astropy.units     as u

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"]
})

def lum_dist(z):
    obj = LambdaCDM(H0=69.6, Om0=0.286, Ode0=1-0.286)
    dist = obj.luminosity_distance(z).value * 3.08567758e24  # cm
    return dist

def lum_flux(f,z):
    return f * 4 * np.pi * lum_dist(z)**2 * 1e-23            # erg s^-1 Hz^-1

def com_dist(z):
    obj = LambdaCDM(H0=69.6, Om0=0.286, Ode0=1-0.286)
    dist = obj.comoving_distance(z).value                    # MPC
    return dist

def com_vol(z_2, z_1=0):
    obj = LambdaCDM(H0=69.6, Om0=0.286, Ode0=1-0.286)
    vol_2 = obj.comoving_volume(z_2).value                   # MPC^3
    vol_1 = obj.comoving_volume(z_1).value                   # MPC^3
    vol = vol_2 - vol_1
    return vol


##########################################################################

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
ID   = []
ID_f = []
Z    = []
SFR  = [[] for i in range(3)]
M    = [[] for i in range(3)]

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

ID   = np.array(ID)
ID_f = np.array(ID_f)
Z    = np.array(Z)
SFR  = np.array(SFR)
M    = np.array(M)


###########################################################################

beautiful_list = ['/data1/astlab09/Beautiful_Aegis.fits',
                  '/data1/astlab09/Beautiful_COSMOS.fits',
                  '/data1/astlab09/Beautiful_GOODSN.fits',
                  '/data1/astlab09/Beautiful_GOODSS.fits',
                  '/data1/astlab09/Beautiful_UDS.fits']

ID_b = []
RA   = []
DEC  = []
    
for beautifulname in beautiful_list:
    beautifulfile = fits.open(beautifulname)
    beauty = beautifulfile[1].data
    beautifulfile.close()
    LB = beauty.shape[0]
    for i in range(LB):
        ID_b.append(beauty[i][0])
        RA.append(beauty[i][2])
        DEC.append(beauty[i][3])
        
ID_b = np.array(ID_b)
RA   = np.array(RA)
DEC  = np.array(DEC)

indexes = []
RA_f    = []
DEC_f   = []
for i in range(len(ID_f)):
    name   = ID_f[i]
    switch = False
    for j in range(len(ID_b)):
        check = ID_b[j]
        if name == check:
            switch = True
            RA_f.append(RA[j])
            DEC_f.append(DEC[j])
    indexes.append(switch)
indexes = np.array(indexes)
ID_f = ID_f[indexes]
ID   = ID[indexes]
Z    = Z[indexes]
new_M   = []
new_SFR = []
for i in range(3):
    new_M.append(M[i][indexes])
    new_SFR.append(SFR[i][indexes])
M   = new_M
SFR = new_SFR


###############################################################


N = len(ID)
f_FUV  = np.zeros(N)
f_NUV  = np.zeros(N)
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


#############################################################


print("cazzo:")
print(len(Z))
print(len(f_FUV[f_FUV != 0.]))


c = 2.99792458e18 # Å / s
lambda_nuv = 2271
lambda_fuv = 1528

fl_FUV = f_FUV * lambda_fuv**2 / c * 10**23
fl_NUV = f_NUV * lambda_nuv**2 / c * 10**23

# ora fl_FUV è in Jy = 10^-23 erg / cm^-2 s^-1 Hz^-1

print("Totale galassie in tutti i fields = {:d}".format(len(Z)))
L_FUV = np.zeros(len(Z))
L_NUV = np.zeros(len(Z))

for i in range(len(Z)):
    L_FUV[i] = lum_flux(fl_FUV[i], Z[i])
    L_NUV[i] = lum_flux(fl_NUV[i], Z[i])

# magnitudine apparente AB = -2.5 log10(flux[Jy] / 3631)    abbiamo cercato come calcolare la magn_app AB in Jy
# plottiamo gli istogrammi per le magnitudini 

m_FUV = np.zeros(len(fl_FUV))
m_NUV = np.zeros(len(fl_NUV))

for i in range(len(fl_FUV)):
    m_FUV[i] = -2.5 * np.log10(fl_FUV[i] / 3631) - 5. * np.log10(lum_dist(Z[i]) / 3.08567758e19)
    m_NUV[i] = -2.5 * np.log10(fl_NUV[i] / 3631) - 5. * np.log10(lum_dist(Z[i]) / 3.08567758e19)

# M_FUV = -2.5 * np.log10(fl_FUV / 3631)
# M_NUV = -2.5 * np.log10(fl_NUV / 3631)

fig , (ax_fuv , ax_nuv) = plt.subplots(2,1,figsize=(15,15))

mag_bins = np.linspace(start=-24,stop=-17,num=11)
ax_fuv.hist(m_FUV, bins = 11, edgecolor='black', color='lightgrey')
ax_nuv.hist(m_NUV, bins = 11, edgecolor='black', color='lightgrey')

# ax_fuv.set_yscale('log')
# ax_nuv.set_yscale('log')
# ax_fuv.set_xlim(-17,-24)
# ax_nuv.set_xlim(-17,-24)
ax_fuv.set_yscale('log')
ax_nuv.set_yscale('log')
ax_fuv.grid(ls=':')  
ax_nuv.grid(ls=':')
ax_fuv.set_xlabel('$m_{FUV}$', fontsize=15)
ax_nuv.set_xlabel('$m_{NUV}$', fontsize=15)
ax_fuv.set_ylabel('$counts$', fontsize=15)
ax_nuv.set_ylabel('$counts$', fontsize=15)


#####################################################################################


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

fig.savefig('/data1/astlab09/0_ALL_FIELDS_histo_luminosity.png', dpi=300)

#######################################################################################

plt.close(fig)

fig , ax = plt.subplots(figsize=(7,7))

# sfr_bibbins = np.linspace(num=7)
sfr_min = np.amin(SFR[1][:])
sfr_max = np.amax(SFR[1][:])
print(sfr_min)
print(sfr_max)
sfr_width = (sfr_max - sfr_min) / 6.
print(sfr_width)

ax.hist(SFR[1][:], bins = 13, edgecolor='black', color='lightgrey')
# ax.set_yscale('log')
ax.set_yscale('log')
ax.set_xlim(None,None)
ax.set_ylim(0,None)
ax.grid(ls=':')
ax.set_xlabel('$log(SFR)$ [$M_{\odot}$ yr$^{-1}$]', fontsize=15)
ax.set_ylabel('$counts$', fontsize=15)

################################################################################


# fissata la galassia: Z[i] , L_FUV[i] , L_NUV[i] , M[:,i] , SFR[:,i] , ID_f[i]
# ordinare per redshift:
    
sort_index = np.argsort(Z)

Z     = Z[sort_index]
L_FUV = L_FUV[sort_index]
L_NUV = L_NUV[sort_index]
m_FUV = m_FUV[sort_index]
m_NUV = m_NUV[sort_index]
ID_f  = ID_f[sort_index]
ID    = ID[sort_index]

for i in range(3):
    M[i]   = M[i][sort_index]
    SFR[i] = SFR[i][sort_index]

    
# Divisione dei dati per bin di redshift
str_Z     = []
str_L_FUV = []
str_L_NUV = []
str_m_FUV = []
str_m_NUV = []
str_ID_f = []
str_ID = []
str_M = [[] for i in range(3)]
str_SFR = [[] for i in range(3)]
str_SSFR = [[] for i in range(3)]

'''
z1 = 1.
z2 = 2.
z3 = 4.
bin_1 = (Z <= z1)                    
bin_2 = (Z > z1) & (Z <= z2)
bin_3 = (Z > z2) & (Z <= z3)
bins  = [bin_1, bin_2, bin_3]
NB = len(bins)
'''
z1 = 0.8
z2 = 1.3
z3 = 1.8
bin_1 = (Z <= z1)                    
bin_2 = (Z > z1) & (Z <= z2)
bin_3 = (Z > z2) & (Z <= z3)
bins  = [bin_1, bin_2, bin_3]
NB = len(bins)

for i in range(len(bins)):
    str_Z.append(Z[bins[i]])
    str_L_FUV.append(L_FUV[bins[i]])
    str_L_NUV.append(L_NUV[bins[i]])
    str_m_FUV.append(m_FUV[bins[i]])
    str_m_NUV.append(m_NUV[bins[i]])
    str_ID_f.append(ID_f[bins[i]])
    str_ID.append(ID[bins[i]])
    for j in range(3):
        str_M[j].append(M[j][bins[i]])
        str_SFR[j].append(SFR[j][bins[i]])
        str_SSFR[j].append(SFR[j][bins[i]] - M[j][bins[i]]) #quantità logaritmiche: - e non /

print(np.amax(m_FUV))
print(np.amax(str_m_FUV[2]))
print(len(str_m_FUV[2]))
print(len(str_SSFR[1][2]))    # stiamo accedendo al percentile 50 del terzo bin di redshift
print(len(str_SFR[1][2]))



##############################################################


d_omega = np.array([3.392190201880897e-05, 
           5.251153102188274e-07, 
           3.559951792354943e-05, 
           1.1785411038838987e-05, 
           1.5116216487179068e-06
          ])

# utilizzando ds9
d_omega = np.array([
           1.044224e-05, 
           1.044603e-05, 
           0.971242e-05, 
           1.362558e-05, 
           1.027001e-05
          ])
d_om = np.sum(d_omega)
# controllare nei paper l'area in °^2 

# volumi dei singoli tronchi di piramide 
V_0 = 1./3. * com_dist(0.3)**3 * d_om
V_1 = 1./3. * com_dist(z1)**3  * d_om - V_0
V_2 = 1./3. * com_dist(z2)**3  * d_om - V_1 - V_0
V_3 = 1./3. * com_dist(z3)**3  * d_om - V_2 - V_1 - V_0

# 1./3. * d_z * ra * d_z * dec * sin(dec) * d_z
# dV = r^2 dr domega
# V = integro in dr = 1/3 r^3 domega

print('V_1:',V_1)
print('V_2:',V_2)
print('V_3:',V_3)

V_1 = com_vol(z_2 = z1, z_1 = 0.3) * d_om / (4. * np.pi) 
V_2 = com_vol(z_2 = z2, z_1 = z1)  * d_om / (4. * np.pi) 
V_3 = com_vol(z_2 = z3, z_1 = z2)  * d_om / (4. * np.pi) 

print('V_1:',V_1)
print('V_2:',V_2)
print('V_3:',V_3)

V = np.array([V_1, V_2, V_3])

# VENGONO BENE!



############################################################################

#####################################################
#####################################################
## LA COSA PIU' SERIA: MAGNITUDE FUNCTION (z BINS) ##
#####################################################
#####################################################

L_sun = 3.826e33 # erg/s

fig , ax = plt.subplots(1, 2, figsize=(12,6), constrained_layout=True)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"]
})

fig.suptitle('Luminosity function - redshift bins', fontsize=15)
ax[0].set_title('Far-UV band', fontsize=15)
ax[1].set_title('Near-UV band', fontsize=15)

# mbins = np.linspace(-24,-12,11)
# centres = np.linspace(-24,-12,21)
# centres = centres[1:20:2] 
mbins = np.array([-17., -17.5, -18., -18.5, -19., -19.5, -20., -20.5, -21., -22., -23., -24.])
mbins = mbins[::-1]
centres = np.zeros(len(mbins) - 1)
for i in range(len(mbins) - 1):
    centres[i] = (mbins[i+1] + mbins[i]) * 0.5
widths  = np.zeros(len(mbins) - 1)
for i in range(len(mbins)-1):
    widths[i] = mbins[i+1] - mbins[i]

density_FUV = []
error_FUV   = []
density_NUV = []
error_NUV   = []
colorstring = ['red','orange','yellow']
labelstring = [   '0.3 $ < z \leq $ {:.1f}'.format(z1)    ,
               '{:.1f} $ < z \leq $ {:.1f}'.format(z1,z2) ,
               '{:.1f} $ < z \leq $ {:.1f}'.format(z2,z3) ] 

for i in range(NB):
    histo_fuv, trash_1 = np.histogram(str_m_FUV[i], bins=mbins)
    histo_nuv, trash_2 = np.histogram(str_m_NUV[i], bins=mbins)
    for j in range(len(histo_fuv)):
        histo_fuv[j] = histo_fuv[j] / widths[j]
        histo_nuv[j] = histo_nuv[j] / widths[j]
    density_FUV.append(histo_fuv / V[i])
    error_FUV.append(np.sqrt(histo_fuv) / V[i])
    density_NUV.append(histo_nuv / V[i])
    error_NUV.append(np.sqrt(histo_nuv) / V[i])
    
    ind_agg     = 0
    fuv_index   = density_FUV[i] != 0.
    nuv_index   = density_NUV[i] != 0.
    fuv_centres = centres[fuv_index]
    nuv_centres = centres[nuv_index]
    fuv_density = density_FUV[i][fuv_index]
    nuv_density = density_NUV[i][nuv_index]
    fuv_error   = error_FUV[i][fuv_index]
    nuv_error   = error_NUV[i][nuv_index]
    fuv_centres = fuv_centres[:np.argmax(fuv_density)+1+ind_agg]
    nuv_centres = nuv_centres[:np.argmax(nuv_density)+1+ind_agg]
    fuv_error   = fuv_error[:np.argmax(fuv_density)+1+ind_agg]
    nuv_error   = nuv_error[:np.argmax(nuv_density)+1+ind_agg]
    fuv_density = fuv_density[:np.argmax(fuv_density)+1+ind_agg]
    nuv_density = nuv_density[:np.argmax(nuv_density)+1+ind_agg]
    
    ax[0].plot(fuv_centres, fuv_density, color=colorstring[i])    
    ax[0].scatter(fuv_centres, fuv_density, color=colorstring[i], label=labelstring[i])
    for j in range(len(fuv_centres)):
        ax[0].vlines(fuv_centres[j], ymin = fuv_density[j] - fuv_error[j], ymax = fuv_density[j] + fuv_error[j], color=colorstring[i])
    ax[0].scatter(fuv_centres, fuv_density - fuv_error, color=colorstring[i], s=50, marker='_')
    ax[0].scatter(fuv_centres, fuv_density + fuv_error, color=colorstring[i], s=50, marker='_')
    
    ax[1].plot(nuv_centres, nuv_density, color=colorstring[i])
    ax[1].scatter(nuv_centres, nuv_density, color=colorstring[i], label=labelstring[i])
    for j in range(len(nuv_centres)):
        ax[1].vlines(nuv_centres[j], ymin = nuv_density[j] - nuv_error[j], ymax = nuv_density[j] + nuv_error[j], color=colorstring[i])
    ax[1].scatter(nuv_centres, nuv_density - nuv_error, color=colorstring[i], s=50, marker='_')
    ax[1].scatter(nuv_centres, nuv_density + nuv_error, color=colorstring[i], s=50, marker='_')
    
    '''
    ax[0].scatter(centres, density_FUV[i], color=colorstring[i], label=labelstring[i])
    for j in range(len(density_NUV[i])):
        ax[0].vlines(centres[j], ymin = density_FUV[i][j] - error_FUV[i][j], ymax = density_FUV[i][j] + error_FUV[i][j], color=colorstring[i])
    ax[0].scatter(centres, density_FUV[i] - error_FUV[i], color=colorstring[i], s=50, marker='_')
    ax[0].scatter(centres, density_FUV[i] + error_FUV[i], color=colorstring[i], s=50, marker='_')
    ax[1].scatter(centres, density_NUV[i], color=colorstring[i], label=labelstring[i])
    for j in range(len(density_NUV[i])):
        ax[1].vlines(centres[j], ymin = density_NUV[i][j] - error_NUV[i][j], ymax = density_NUV[i][j] + error_NUV[i][j], color=colorstring[i])
    ax[1].scatter(centres, density_NUV[i] - error_NUV[i], color=colorstring[i], s=50, marker='_')
    ax[1].scatter(centres, density_NUV[i] + error_NUV[i], color=colorstring[i], s=50, marker='_')
    '''

for j in range(2):
    ax[j].set_yscale('log')
    ax[j].set_xlim(-17,-24)
    # ax[j].set_ylim(1e-4,1e3)
    ax[j].set_ylim(1e-7,1e-2)
    ax[j].grid(ls=':',which='both')
    ax[j].set_ylabel('$n$ [Mpc$^{-3}$ mag$^{-1}$]', fontsize=15) 
    ax[j].legend(frameon=True)
ax[0].set_xlabel('$M_{FUV}$ [mag]', fontsize=15)
ax[1].set_xlabel('$M_{NUV}$ [mag]', fontsize=15)

if (z1 == 1.) and (z2 == 2.) and (z3 == 4.):
    fig.savefig('/data1/astlab09/ALL_FIELDS_magnitude_function_z_bins_large.png', dpi=300)
    fig.savefig('/data1/astlab09/ALL_FIELDS_magnitude_function_z_bins_large.eps')
    fig.savefig('/data1/astlab09/ALL_FIELDS_magnitude_function_z_bins_large.pdf')
elif (z1 == 0.8) and (z2 == 1.3) and (z3 == 1.8):
    fig.savefig('/data1/astlab09/ALL_FIELDS_magnitude_function_z_bins_tight.png', dpi=300)
    fig.savefig('/data1/astlab09/ALL_FIELDS_magnitude_function_z_bins_tight.eps')
    fig.savefig('/data1/astlab09/ALL_FIELDS_magnitude_function_z_bins_tight.pdf')



#################################################################################


######################################################
######################################################
## LA COSA PIU' SERIA: LUMINOSITY FUNCTION (z BINS) ##
######################################################
######################################################

L_sun = 3.826e33 # erg/s

fig , ax = plt.subplots(1, 2, figsize=(12,6), constrained_layout=True)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"]
})

fig.suptitle('Luminosity function - redshift bins', fontsize=15)
ax[0].set_title('Far-UV band', fontsize=15)
ax[1].set_title('Near-UV band', fontsize=15)

logbins = np.logspace(25,30,11) / L_sun
centres = np.logspace(25,30,21) / L_sun
centres = centres[1:20:2] 
widths  = np.zeros(10)
for i in range(len(logbins)-1):
    widths[i] = logbins[i+1] - logbins[i]

density_FUV = []
density_NUV  = []
colorstring = ['red','orange','yellow']
labelstring = [   '0.3 $ < z \leq $ {:.1f}'.format(z1)    ,
               '{:.1f} $ < z \leq $ {:.1f}'.format(z1,z2) ,
               '{:.1f} $ < z \leq $ {:.1f}'.format(z2,z3) ] 

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
    ax[j].legend(frameon=True)
ax[0].set_xlabel('$L_{FUV}$ [$L_{\odot}$ Hz$^{-1}$]', fontsize=15)
ax[1].set_xlabel('$L_{NUV}$ [$L_{\odot}$ Hz$^{-1}$]', fontsize=15)


if (z1 == 1.) and (z2 == 2.) and (z3 == 4.):
    fig.savefig('/data1/astlab09/ALL_FIELDS_luminosity_function_z_bins_large.png', dpi=300)
    fig.savefig('/data1/astlab09/ALL_FIELDS_luminosity_function_z_bins_large.eps')
    fig.savefig('/data1/astlab09/ALL_FIELDS_luminosity_function_z_bins_large.pdf')
elif (z1 == 0.8) and (z2 == 1.3) and (z3 == 1.8):
    fig.savefig('/data1/astlab09/ALL_FIELDS_luminosity_function_z_bins_tight.png', dpi=300)
    fig.savefig('/data1/astlab09/ALL_FIELDS_luminosity_function_z_bins_tight.eps')
    fig.savefig('/data1/astlab09/ALL_FIELDS_luminosity_function_z_bins_tight.pdf')



#############################################################################


###############################################
###############################################
## LA COSA PIU' SERIA: SFR FUNCTION (z BINS) ##
###############################################
###############################################

fig , ax = plt.subplots(figsize=(6,6), constrained_layout=True)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"]
})

ax.set_title('Star formation rate function - redshift bins', fontsize=15)

mbins = np.linspace(sfr_min,sfr_max,21)
# centres = np.linspace(sfr_min,sfr_max,15)
# centres = centres[1:14:2] 
# mbins = np.array([-17., -17.5, -18., -18.5, -19., -19.5, -20., -20.5, -21., -22., -23., -24.])
# mbins = mbins[::-1]
centres = np.zeros(len(mbins) - 1)
for i in range(len(mbins) - 1):
    centres[i] = (mbins[i+1] + mbins[i]) * 0.5
widths  = np.zeros(len(mbins) - 1)
for i in range(len(mbins)-1):
    widths[i] = mbins[i+1] - mbins[i]

density = []
error   = []
colorstring = ['midnightblue','dodgerblue','lightskyblue']
labelstring = [   '0.3 $ < z \leq $ {:.1f}'.format(z1)    ,
               '{:.1f} $ < z \leq $ {:.1f}'.format(z1,z2) ,
               '{:.1f} $ < z \leq $ {:.1f}'.format(z2,z3) ] 

for i in range(3):
    histo, trash = np.histogram(str_SFR[1][i], bins=mbins)
    for j in range(len(histo)):
        histo[j] = histo[j] / widths[j]
    density.append(histo / V[i])
    error.append(np.sqrt(histo) / V[i])
    
    s_index   = density[i] != 0.
    s_centres = centres[s_index]
    s_density = density[i][s_index]
    s_error   = error[i][s_index]
    
    s_centres = s_centres[np.argmax(s_density):]
    s_error   = s_error[np.argmax(s_density):]
    s_density = s_density[np.argmax(s_density):]
    
    ax.plot(s_centres, s_density, color=colorstring[i])    
    ax.scatter(s_centres, s_density, color=colorstring[i], label=labelstring[i])
    for j in range(len(s_centres)):
        ax.vlines(s_centres[j], ymin = s_density[j] - s_error[j], ymax = s_density[j] + s_error[j], color=colorstring[i])
    ax.scatter(s_centres, s_density - s_error, color=colorstring[i], s=50, marker='_')
    ax.scatter(s_centres, s_density + s_error, color=colorstring[i], s=50, marker='_')
    
ax.set_yscale('log')
ax.set_xlim(-1,sfr_max)
# ax[j].set_ylim(1e-4,1e3)
# ax.set_ylim(2e-7,2e-3)
ax.grid(ls=':',which='both')
ax.set_ylabel('$n$ [Mpc$^{-3}$ ($M_{\odot}$ yr$^{-1}$)$^{-1}$]', fontsize=15) 
ax.legend(frameon=True)
ax.set_xlabel('$log(SFR)$ [$M_{\odot}$ yr$^{-1}$]', fontsize=15)

if (z1 == 1.) and (z2 == 2.) and (z3 == 4.):
    fig.savefig('/data1/astlab09/ALL_FIELDS_SFR_function_z_bins_large.png', dpi=300)
    fig.savefig('/data1/astlab09/ALL_FIELDS_SFR_function_z_bins_large.eps')
    fig.savefig('/data1/astlab09/ALL_FIELDS_SFR_function_z_bins_large.pdf')
elif (z1 == 0.8) and (z2 == 1.3) and (z3 == 1.8):
    fig.savefig('/data1/astlab09/ALL_FIELDS_SFR_function_z_bins_tight.png', dpi=300)
    fig.savefig('/data1/astlab09/ALL_FIELDS_SFR_function_z_bins_tight.eps')
    fig.savefig('/data1/astlab09/ALL_FIELDS_SFR_function_z_bins_tight.pdf')


##############################################################################

################################################
################################################
## LA COSA PIU' SERIA: SSFR FUNCTION (z BINS) ##
################################################
################################################

ssfr_min = np.amin(SFR[1][:] - M[1][:])
ssfr_max = np.amax(SFR[1][:] - M[1][:])
print(ssfr_min)
print(ssfr_max)
ssfr_width = (ssfr_max - ssfr_min) / 6.
print(ssfr_width)


fig , ax = plt.subplots(figsize=(6,6), constrained_layout=True)

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Times New Roman"]
})

ax.set_title('Specific star formation rate function - redshift bins', fontsize=15)

mbins = np.linspace(ssfr_min,ssfr_max,21)
# centres = np.linspace(sfr_min,sfr_max,15)
# centres = centres[1:14:2] 
# mbins = np.array([-17., -17.5, -18., -18.5, -19., -19.5, -20., -20.5, -21., -22., -23., -24.])
# mbins = mbins[::-1]
centres = np.zeros(len(mbins) - 1)
for i in range(len(mbins) - 1):
    centres[i] = (mbins[i+1] + mbins[i]) * 0.5
widths  = np.zeros(len(mbins) - 1)
for i in range(len(mbins)-1):
    widths[i] = mbins[i+1] - mbins[i]

density = []
error   = []
colorstring = ['midnightblue','dodgerblue','lightskyblue']
labelstring = [   '0.3 $ < z \leq $ {:.1f}'.format(z1)    ,
               '{:.1f} $ < z \leq $ {:.1f}'.format(z1,z2) ,
               '{:.1f} $ < z \leq $ {:.1f}'.format(z2,z3) ] 

for i in range(3):
    histo, trash = np.histogram(str_SSFR[1][i], bins=mbins)
    for j in range(len(histo)):
        histo[j] = histo[j] / widths[j]
    density.append(histo / V[i])
    error.append(np.sqrt(histo) / V[i])
    
    s_index   = density[i] != 0.
    s_centres = centres[s_index]
    s_density = density[i][s_index]
    s_error   = error[i][s_index]
    '''
    s_centres = s_centres[np.argmax(s_density):]
    s_error   = s_error[np.argmax(s_density):]
    s_density = s_density[np.argmax(s_density):]
    '''
    ax.plot(s_centres, s_density, color=colorstring[i])    
    ax.scatter(s_centres, s_density, color=colorstring[i], label=labelstring[i])
    for j in range(len(s_centres)):
        ax.vlines(s_centres[j], ymin = s_density[j] - s_error[j], ymax = s_density[j] + s_error[j], color=colorstring[i])
    ax.scatter(s_centres, s_density - s_error, color=colorstring[i], s=50, marker='_')
    ax.scatter(s_centres, s_density + s_error, color=colorstring[i], s=50, marker='_')
    
ax.set_yscale('log')
ax.set_xlim(-14,ssfr_max) # ssfr_min
# ax[j].set_ylim(1e-4,1e3)
# ax.set_ylim(2e-7,2e-3)
ax.grid(ls=':',which='both')
ax.set_ylabel('$n$ [Mpc$^{-3}$ ($M_{\odot}$ yr$^{-1}$)$^{-1}$]', fontsize=15) 
ax.legend(frameon=True)
ax.set_xlabel('$log(SSFR)$ [yr$^{-1}$]', fontsize=15)

if (z1 == 1.) and (z2 == 2.) and (z3 == 4.):
    fig.savefig('/data1/astlab09/ALL_FIELDS_SSFR_function_z_bins_large.png', dpi=300)
    fig.savefig('/data1/astlab09/ALL_FIELDS_SSFR_function_z_bins_large.eps')
    fig.savefig('/data1/astlab09/ALL_FIELDS_SSFR_function_z_bins_large.pdf')
elif (z1 == 0.8) and (z2 == 1.3) and (z3 == 1.8):
    fig.savefig('/data1/astlab09/ALL_FIELDS_SSFR_function_z_bins_tight.png', dpi=300)
    fig.savefig('/data1/astlab09/ALL_FIELDS_SSFR_function_z_bins_tight.eps')
    fig.savefig('/data1/astlab09/ALL_FIELDS_SSFR_function_z_bins_tight.pdf')


