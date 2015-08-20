#!/usr/bin/env python                                                                                                                                                                 
###################                                                                                                                                                               
# Shany Danieli   #
# Yale University #                                                                                                                                                               
# 9/23/2014       #
###################

# 1/23/2015 - adding the calculation of the covariance matrix for the errors
                                                                                                                                                                   
##############################################################################                                                                                                    
# The program calculates and plots the histogram of the luminosity           #
# function per unit volume.                                                  #
# The mock was taken from Andrew Hearin and produced using the Bolshoi       #
# N-body simulation.                                                         #
#                                                                            #                                                                                                    
# Required input:                                                            #
# To run type: python luminosityFun.py                                       #                                                                                                    
#                                                                            #
##############################################################################   

# Import python modules
from pylab import *
from numpy import *
from math import *
import matplotlib.pyplot as plt
from scipy import interpolate
import h5py



V = 250.0**3


"""
Data processing - Galaxy Luminosity Function
"""

# Read in the data from the file "conditional_abundance_matching_Mr_gr_model.dat"
#data = loadtxt("../../Data/conditional_abundance_matching_Mr_gr_model.dat",float) #abundance matching mock with assembly bias
data = loadtxt("../../Data/erased_assembly_bias_Mr_gr_model.dat",float) #abundance matching mock with erased assembly bias
#data = loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/erased_assembly_bias_Mr_gr_model_downsized.dat",float) #abundance matching mock with erased assembly bias, downsized to Mhost=[11.75-14.25]

#data2 = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/sham_mock.dat",float) # with assembly bias
#data2 = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock_fixed.dat",float) # no assembly bias
data2 = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock.dat",float) # no assembly bias
#data2 = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock_downsized.dat",float) # no assembly bias, downsized to Mhost=[11.75-14.25]

#halo_ID = data2[:,0] # Halo ID
#UPID = data2[:,1] # if UPID == -1 -> central, if UPID > 0 -> satellite



"""
Mhost = data2[:,2] # Host mass in units of (Msun/h)
absMag = data[:,7]
pos = data[:,range(1,4)] # the position of each galaxy in the box
N = len(absMag)
"""

"""
# Downsizing the sample
Mhost_tmp = data2[:,2] # Host mass in units of (Msun/h)
Mhost = np.asarray([Mhost_tmp[i] for i in range(len(Mhost_tmp)) if (11.75-dm/2)<log10(Mhost_tmp[i])<=(14.25+dm/2)])
absMag_tmp = data[:,7]
absMag = np.asarray([absMag_tmp[i] for i in range(len(absMag_tmp)) if (11.75-dm/2)<log10(Mhost_tmp[i])<=(14.25+dm/2)])
pos_tmp = data[:,range(1,4)] # the position of each galaxy in the box
pos = np.asarray([pos_tmp[i] for i in range(len(pos_tmp)) if (11.75-dm/2)<log10(Mhost_tmp[i])<=(14.25+dm/2)])
upid_tmp = data2[:,1] # Host mass in units of (Msun/h)
upid = np.asarray([upid_tmp[i] for i in range(len(upid_tmp)) if (11.75-dm/2)<log10(Mhost_tmp[i])<=(14.25+dm/2)])
N = len(absMag)


# Split sample into central and satellite galaxies
mag_cen = []
mag_sat = []
for i in range(len(absMag)):
    if upid[i]<0: # central galaxies
        mag_cen.append(absMag[i])
    else: # satellite galaxies
        mag_sat.append(absMag[i])

"""
# full sample
#data1_cen = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/splitted/cen_mag.txt",float)
#data2_cen = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/splitted/cen.txt",float)

# downsized sample
data1_cen = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/downsized/cen_mag_ds.txt",float)
data2_cen = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/downsized/cen_ds.txt",float)


abs_mag_cen = data1_cen[:,7]
Msun = 4.76
L_data_cen = ((Msun-abs_mag_cen)/2.5) # L is log(L/L_sun)
host_mass_cen = data2_cen[:,2]
upid_cen = data2_cen[:,0]

#only satellites
# full sample
#data1_sat = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/splitted/sat_mag.txt",float)
#data2_sat = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/splitted/sat.txt",float)

# downsized sample
data1_sat = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/downsized/sat_mag_ds.txt",float)
data2_sat = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/downsized/sat_ds.txt",float)


abs_mag_sat = data1_sat[:,7]
L_data_sat = ((Msun-abs_mag_sat)/2.5) # L is log(L/L_sun)
host_mass_sat = data2_sat[:,2]
upid_sat = data2_sat[:,1]


"""
Luminosity function calculation
"""


# Luminosity bins
L_min = 9.5 # 9.8	
L_max = 10.7	
#L_min = 9.7
#L_max = 11.2	
dL = 0.15
L_bins = np.linspace(L_min,L_max,(L_max-L_min)/dL+1)
L_centers = [x+dL/2 for x in L_bins]
del L_centers[-1]

# satellite galaxies
(counts, L_bins, patches) = hist(L_data_sat,L_bins) 
dndlogL_sat = counts/(dL*V)
plt.close()

# central galaxies
(counts, L_bins, patches) = hist(L_data_cen,L_bins) 
dndlogL_cen = counts/(dL*V)
plt.close()



"""
# Calculation of the luminosity function
Vbox = 250.0**3 # The volume of the box
delta_M = 0.15 # size of the absolute magnitude bins
bins = [x*0.15-22 for x in range(0,21)]
(n, bins, patches) = hist(absMag,bins,color='w') # 0.15 bins from -19 to -22
Nbins = len(bins)
# Absolute magnitude bins
magnitude = [round(x-delta_M/2,3) for x in bins]
del magnitude[0]
Msun = 4.76
L_data = ((Msun-np.asarray(absMag))/2.5) # L is log(L/L_sun)

dndMag = [x/(Vbox*delta_M) for x in n] # the luminosity function as a function of absolute magnitude
# Luminosity bins
#L_min = 9.5 # 9.8	
#L_max = 11.2	

L_min = 9.5 
L_max = 10.9	

dL = 0.15

L_bins = np.linspace(L_min,L_max,(L_max-L_min)/dL+1)
L_centers = [x+dL/2 for x in L_bins]
del L_centers[-1]
# Calculates the luminosity function as a function of the luminosity
(k, L_bins, patches) = hist(L_data,L_bins,color='w') 
plt.close()
dndlogL = [x/(Vbox*dL) for x in k] # the luminosity function as a function of luminosity



# Luminosity function for centrals and satellites
L_cen = ((Msun-np.asarray(mag_cen))/2.5) # L is log(L/L_sun)
L_sat = ((Msun-np.asarray(mag_sat))/2.5) # L is log(L/L_sun)
(counts_cen, L_bins, patches1) = hist(L_cen,L_bins,color='w')
plt.close()
dndlogL_cen = np.asarray([x/(Vbox*dL) for x in counts_cen])
(counts_sat, L_bins, patches2) = hist(L_sat,L_bins,color='w')
plt.close()
dndlogL_sat = np.asarray([x/(Vbox*dL) for x in counts_sat])
"""
dndlogL = dndlogL_cen + dndlogL_sat

print "Computed directly"
print "dndlogL_cen:"
print dndlogL_cen
print "\n"
print "dndlogL_sat:"
print dndlogL_sat
print "\n\n"


"""
Jackknife errors calculation
"""

"""
# Calculates the errors using the Jackknife method
Lbox = 250.0
V = 250.0**3
Ndivs = 5
dl = Lbox/Ndivs
n_subBox = Ndivs*Ndivs*Ndivs # The number of sub volumes
V_subBox = Vbox - Vbox/n_subBox


# Indices for the galaxies positions
index = asarray([floor(pos[i,0]/dl) + (floor(pos[i,1]/dl)*Ndivs) + (floor(pos[i,2]/dl)*Ndivs*Ndivs) + 1 for i in range(N)]) # index for the position of particle2
subBox = []

#absMag_all = []  # 2D list for the values of the luminosity function for each one of the bins and for each sub volume
#absMag_all.append(dndMag)

L_data_all = []
L_data_all.append(dndlogL)

for k in range(1,n_subBox+1): # run over the sub-samples
        for i in range(0,N): # runs over all the points (galaxies)
                if (index[i] != k): # the point is inside the sub-box
#                    subBox.append(absMag[i]) # then add to sub-box list
                    subBox.append(L_data[i]) # then add to sub-box list
#        (n_subVol, bins_subBox, patches_subBox) = hist(subBox, bins)
        (n_subVol, bins_subBox, patches_subBox) = hist(subBox, L_bins)
#        dndMag_sub = [x/(V_subBox*delta_M) for x in n_subVol] # the luminosity function
        dndlogL_sub = [x/(V_subBox*dL) for x in n_subVol] # the luminosity function
#        absMag_all.append(dndMag_sub)
        L_data_all.append(dndlogL_sub)
        subBox = []





full = L_data_all[0] # the luminosity function for the full sample
full = asarray(full)
subSamples = L_data_all[1:] # the luminosity function for the sub-samples
subSamples = asarray(subSamples)
after_subtraction = subSamples - np.mean(subSamples,axis=0)
n_subBox = float(n_subBox) # The number of sub volumes
"""


"""
full = absMag_all[0] # the luminosity function for the full sample
full = asarray(full)
subSamples = absMag_all[1:] # the luminosity function for the sub-samples
"""


"""
Nr = full.shape[0] # Nr is the number of radial bins
cov = np.zeros((Nr,Nr)) # 2D array that keeps the covariance matrix
savetxt("Output/after_subtraction_lumi.dat", after_subtraction,fmt='%10.11f') # saving the after_subtraction matrix for the calculation of the full covariance matrix

tmp = 0
for i in range(Nr):
    for j in range(Nr):
        for k in range(int(n_subBox)):
            tmp = tmp + after_subtraction[k,i]*after_subtraction[k,j]
        cov[i,j] = (((n_subBox-1)/n_subBox)*tmp)**0.5
        tmp = 0
#savetxt("Output/lumi_cov.dat", cov,fmt='%10.9f')
cov_reversed = np.zeros((Nr,Nr)) # 2D array that keeps the covariance matrix
for p in range(0,Nr):
    for t in range(0,Nr):
        cov_reversed[p,t] = cov[Nr-1-p][Nr-1-t] # 2D array that keeps the covariance matrix
savetxt("Output/lumi_cov.dat", cov_reversed,fmt='%10.9f')
"""

"""
squared = after_subtraction**2
error2 = ((n_subBox-1)/n_subBox)*squared.sum(axis=0)
errors = error2**0.5
# reverse the order
#magnitude.reverse()
#dndMag.reverse()
#dndlogL.reverse()
errors = errors.tolist()
#errors.reverse()
plt.close() 
"""




"""
for f1,f2,f3 in zip(magnitude,dndMag,errors):
    print (str(f1)+ "     " + str(log(f2)) + "     " + str(f3))
    print ("\n")
"""

"""
# convert the absolute magnitude to luminosity in order to compare later 
Msun = 4.76 # The Sun's absolute magnitude
L = ((Msun-np.asarray(magnitude))/2.5) # L is log(L/L_sun)
L_not_log = [10**x for x in L] # this is L/L_sun
lumi_log = [log(n,10) for n in dndMag]
"""

"""
# Write into file
f = open("Output/direct_calc.dat","w")
for f1,f2,f3 in zip(L,dndMag,errors):
    f.write(str(f1)+ "     " + str(f2) + "     " + str(f3))
    f.write("\n")
f.close()
"""


# importing the CLF computed in a different file
clf_cen = loadtxt("/Users/si224/Documents/2014/vanDenBosch/Code/CLF/output/clf_cen.dat",float)[:,:]
clf_sat = loadtxt("/Users/si224/Documents/2014/vanDenBosch/Code/CLF/output/clf_sat.dat",float)[:,:]



# Host mass bins
num_bins = 15
Mmin = 11.75 # left edge
Mmax = 14.25 # right edge
dlogM = (Mmax-Mmin)/num_bins
Mass_bins = np.arange(Mmin,Mmax+dlogM/2,dlogM)
Mass_centers = (Mass_bins+dlogM/2).tolist()
del Mass_centers[-1]


"""
# Host mass bins
num_bins = 20
#Mmin = 10.4
#Mmax = 15.0
Mmin = 11.85
Mmax = 14.15
dlogM = (Mmax-Mmin)/(num_bins-1)
Mass_centers = np.arange(Mmin,Mmax+dlogM/2,dlogM)
Mass_bins = (Mass_centers-dlogM/2).tolist()
Mass_bins.append(Mass_bins[-1]+dlogM)
"""

# Analytical computation of the luminosity function
# Uses: lumi_fun = sum_over_M{CLF(L|M)*n(M)}


# Analytical halo mass function
masses_analytic = loadtxt("/Users/si224/Documents/2014/vanDenBosch/Code/mass_function/MF_Code_Tinker/tinker_"+str(num_bins)+"_bins.dndM",float)[:,0]
hmf_analytic = loadtxt("/Users/si224/Documents/2014/vanDenBosch/Code/mass_function/MF_Code_Tinker/tinker_"+str(num_bins)+"_bins.dndM",float)[:,1]
factor = log(10)*masses_analytic  # Converting from dndM to dndlogM
hmf_analytic = factor*hmf_analytic
dlogM_model = log10(masses_analytic[1])-log10(masses_analytic[0])

# Analytical 
lumi_clf_before_sum = [clf_cen[i]*hmf_analytic[i] for i in range(len(hmf_analytic))] + [clf_sat[i]*hmf_analytic[i] for i in range(len(hmf_analytic))]
lumi_clf = sum(lumi_clf_before_sum,axis=0)*dlogM_model

lumi_clf_cen = sum([clf_cen[i]*hmf_analytic[i] for i in range(len(hmf_analytic))],axis=0)*dlogM_model
lumi_clf_sat = sum([clf_sat[i]*hmf_analytic[i] for i in range(len(hmf_analytic))],axis=0)*dlogM_model




"""
# Empirical halo mass function computed from the simulation
hmf_data =loadtxt("/Users/si224/Documents/2014/vanDenBosch/Code/mass_function/hmf_data.dat",float)
#dlogM_data = (Mmax-Mmin)/(num_bins-1)
dlogM_data = (Mmax-Mmin)/num_bins
#bin_centers_data = np.arange(Mmin,Mmax+dlogM_data/2,dlogM_data)


# Empirical
lumi_clf_before_sum = [clf_cen[i]*hmf_data[i] for i in range(len(hmf_data))] + [clf_sat[i]*hmf_data[i] for i in range(len(hmf_data))]
lumi_clf = sum(lumi_clf_before_sum,axis=0)*dlogM_data

lumi_clf_cen = sum([clf_cen[i]*hmf_data[i] for i in range(len(hmf_data))],axis=0)*dlogM_data
lumi_clf_sat = sum([clf_sat[i]*hmf_data[i] for i in range(len(hmf_data))],axis=0)*dlogM_data
"""




diff = (np.log10(dndlogL)-np.log10(lumi_clf))/np.log10(dndlogL)
#diff = (np.log10(dndlogL)-np.log10(lumi_clf))/np.log10(lumi_clf)
print "diff:"
print 100*diff
diff = (dndlogL-lumi_clf)/dndlogL
print "diff:"
print 100*diff

#diff_cen = (dndlogL_cen-lumi_clf_cen)/dndlogL_cen
#diff_sat = (dndlogL_sat-lumi_clf_sat)/dndlogL_sat

diff_cen = (np.log10(dndlogL_cen)-np.log10(lumi_clf_cen))/np.log10(dndlogL_cen)
diff_sat = (np.log10(dndlogL_sat)-np.log10(lumi_clf_sat))/np.log10(dndlogL_sat)
#diff_sat = np.delete(diff_sat, len(diff_sat)-1)


"""
lumi_clf_sat = lumi_clf_sat.tolist()
lumi_clf_cen = lumi_clf_cen.tolist()

del lumi_clf_sat[-1]
del dndlogL_sat[-1]
"""




plt.figure
plt.subplot(2,1,1)
plt.scatter(L_centers,[log(x,10) for x in lumi_clf_cen],color = 'b',marker='D',label='Using CLF')
plt.scatter(L_centers,[log(x,10) for x in dndlogL_cen],color = 'r',marker='o',label='Directly')
plt.legend()
plt.xlabel(r'L[L$_{\odot}$h$^{-2}$]',fontsize=15)
plt.ylabel(r'log(L$\Phi$(L))[h/Mpc]$^{3}$',fontsize=15)
plt.title('central galaxies')
plt.subplot(2,1,2)
plt.plot(L_centers,diff_cen,linestyle='None',marker='o',markersize=5,markeredgecolor='k', markerfacecolor='m')
plt.xlabel(r'L[L$_{\odot}$h$^{-2}$]',fontsize=15)
plt.ylabel(r'(Simulation-Model)/Model',fontsize=15)
plt.ylim(-0.4,0.4)
plt.axhline(linewidth=0.5, color='k')
plt.axhspan(-0.05, 0.05, facecolor='0.5',linewidth=0)
plt.axhspan(-0.05, -0.1, facecolor='0.7',linewidth=0)
plt.axhspan(0.05, 0.1, facecolor='0.7',linewidth=0)
plt.axhspan(0.1, 0.2, facecolor='0.9',linewidth=0)
plt.axhspan(-0.1, -0.2, facecolor='0.9',linewidth=0)
plt.show()
plt.close()



plt.figure
plt.subplot(2,1,1)
L_centers_sat = [x for x in L_centers]
#del L_centers_sat[-1]
plt.scatter(L_centers_sat,[log(x,10) for x in lumi_clf_sat],color = 'g',marker='o',label='Using CLF')
plt.scatter(L_centers_sat,[log(x,10) for x in dndlogL_sat],color = 'm',marker='o',label='Directly')
plt.legend()
plt.xlabel(r'L[L$_{\odot}$h$^{-2}$]',fontsize=15)
plt.ylabel(r'log(L$\Phi$(L))[h/Mpc]$^{3}$',fontsize=15)
plt.title('satellite galaxies')
plt.subplot(2,1,2)
plt.plot(L_centers_sat,diff_sat,linestyle='None',marker='o',markersize=5,markeredgecolor='k', markerfacecolor='m')
plt.xlabel(r'L[L$_{\odot}$h$^{-2}$]',fontsize=15)
plt.ylabel(r'[$\Phi_{direct}$-$\Phi_{CLF}$]/$\Phi_{direct}$',fontsize=15)
plt.ylim(-0.4,0.4)
plt.axhline(linewidth=0.5, color='k')
plt.axhspan(-0.05, 0.05, facecolor='0.5',linewidth=0)
plt.axhspan(-0.05, -0.1, facecolor='0.7',linewidth=0)
plt.axhspan(0.05, 0.1, facecolor='0.7',linewidth=0)
plt.axhspan(0.1, 0.2, facecolor='0.9',linewidth=0)
plt.axhspan(-0.1, -0.2, facecolor='0.9',linewidth=0)
plt.legend()
plt.show()




"""
Interpolating the CLF
"""

"""
f = interpolate.interp1d(L_centers, lumi_clf,'cubic')
L_centers_new = np.arange(L_centers[0],L_centers[-1],0.001)
lumi_clf_new = f(L_centers_new)
"""


"""
Plotting
"""

# plot the luminosity function
# with the analytical calculation from the CLF and hmf
#errors_log = [x/y for x,y in zip(errors,dndlogL)]
#plt.figure()
plt.subplot(2,1,1)
#plt.errorbar(L_centers,[log(x,10) for x in dndlogL], yerr=errors_log, ls='none',color = 'b', elinewidth=2,marker='o',label='Simulation')
plt.plot(L_centers,[log(x,10) for x in dndlogL],linestyle='none',color = 'b',marker='D', markersize=7.0,markeredgewidth=0.5,label='Direct',zorder = 1)
plt.plot(L_centers,[log(x,10) for x in lumi_clf],linestyle='none',color='c',markeredgecolor='k',markeredgewidth=0.5,markersize=5.0,marker='o',label='Using CLF',zorder = 2)
plt.xlabel(r'L[L$_{\odot}$h$^{-2}$]',fontsize=15)
plt.ylabel(r'log(L$\Phi$(L))[h/Mpc]$^{3}$',fontsize=15)
plt.legend()
plt.subplot(2,1,2)
plt.plot(L_centers,diff,linestyle='None',marker='o',markersize=5,markeredgecolor='k', markerfacecolor='m')
plt.xlabel(r'L[L$_{\odot}$h$^{-2}$]',fontsize=15)
plt.ylabel(r'[$\Phi_{direct}$-$\Phi_{CLF}$]/$\Phi_{direct}$',fontsize=15)
plt.ylim(-0.4,0.4)
plt.axhline(linewidth=0.5, color='k')
plt.axhspan(-0.05, 0.05, facecolor='0.5',linewidth=0)
plt.axhspan(-0.05, -0.1, facecolor='0.7',linewidth=0)
plt.axhspan(0.05, 0.1, facecolor='0.7',linewidth=0)
plt.axhspan(0.1, 0.2, facecolor='0.9',linewidth=0)
plt.axhspan(-0.1, -0.2, facecolor='0.9',linewidth=0)
plt.legend()
plt.savefig('directVsCLF.png')
plt.show()
plt.close()


#plt.errorbar(magnitude,lumi_log, yerr=errors_log, capsize=4, ls='none',color = 'k', elinewidth=1,marker='x',label='Direct from simulation')
#plt.scatter(magnitude,[log(x,10) for x in lumi_model],color = 'b',marker='o',label='Model(Tinker)')
#plt.plot(L_centers_new,[log(x,10) for x in lumi_clf_new],'-',c='k')
#errors_log = [x/y for x,y in zip(errors,dndMag)]
#plt.xlabel(r'$M_r-5\logh$',fontsize=15)
#plt.ylabel(r'$\log(\Phi[h^{3}Mpc^{-3}mag^{-1}])$',fontsize=15)
#plt.gca().invert_xaxis()
#plt.savefig('directVsCLF_11.75-14.25_full_sample.png')
#plt.gca().invert_xaxis()










####################################################################################
# Code that I commented along the way.											   #
#																				   #
#																				   #
#																				   #
#																				   #
####################################################################################

"""
Data processing for the stellar mass function (SMF)
"""
"""
Vbox = 15625000 # The volume of the box
f = h5py.File('/Users/si224/Documents/2014/vanDenBosch/Data/Duncan/sm_9.5_s0.0_sfr_c-0.25_250.hdf5', 'r') #open mock in 'read' mode
catalogue = 'sm_9.5_s0.0_sfr_c-0.25_250'
mock = f.get(catalogue)

halo_ID = mock['id'] 
UPID = mock['pid']
Mhost = mock['mvir'] 

print "min and max:"
print max(Mhost)
print min(Mhost)
L_data = mock['Mstar'] 
pos = [mock['x'],mock['y'],mock['z']]
N = len(halo_ID)
"""

"""
stellar mass function calculation
"""
"""
# Calculation of the stellar mass function
# Stellar mass bins
L_min = 9.5
L_max = 11.6	
dL = 0.15
L_bins = np.linspace(L_min,L_max,(L_max-L_min)/dL+1)
L_centers = [x+dL/2 for x in L_bins]
del L_centers[-1]
# Calculates the luminosity function as a function of the luminosity
(k, L_bins, patches) = hist(L_data,L_bins,color='w') 
plt.close()
dndlogL = [x/(Vbox*dL) for x in k] # the luminosity function as a function of luminosity


L_cen = []
L_sat = []
for i in range(len(L_data)):
    if UPID[i]<0:
        L_cen.append(L_data[i])
    else:
        L_sat.append(L_data[i])

# Luminosity function for centrals and satellites
(counts_cen, L_bins, patches) = hist(L_cen,L_bins,color='w')
plt.close()
dndlogL_cen = [x/(Vbox*dL) for x in counts_cen]
(counts_sat, L_bins, patches) = hist(L_sat,L_bins,color='w')
plt.close()
dndlogL_sat = [x/(Vbox*dL) for x in counts_sat]

print "dndlogL:"
print np.log10(dndlogL)




#plt.scatter(L_centers,[log(x,10) for x in dndlogL])#, ls='none',color = 'b', elinewidth=2,marker='o',label='stellar mass function - simulation')
#plt.show()
"""





"""
clf_cen_transpose = asarray(clf_cen).T.tolist()
savetxt("output/clf_cen_transpose.dat", clf_cen_transpose)
clf_sat_transpose = asarray(clf_sat).T.tolist()
savetxt("output/clf_sat_transpose.dat", clf_sat_transpose)
"""

"""
lumi_cen = sum(clf_cen_transpose*hmf_analytic, axis=1) # need to multiply by dlogM (the bin width)
lumi_sat = sum(clf_sat_transpose*hmf_analytic, axis=1)
lumi_model = 0.4*(lumi_cen + lumi_sat)*dlogM
#lumi = 0.4*(lumi_cen + lumi_sat)
diff = (dndMag-lumi_model)/dndMag
"""

""" commented on 7/16/2015
lumi_log = [log(n,10) for n in dndMag]
errors_log = [x/y for x,y in zip(errors,dndMag)]
plt.errorbar(magnitude,lumi_log, yerr=errors_log, capsize=4, ls='none',color = 'k', elinewidth=1,marker='x')
plt.xlim(-19,-22)
plt.xlabel(r'$M_r-5\logh$',fontsize=15)
plt.ylabel(r'$\log(\Phi/[h^3Mpc^-3mag^-1])$',fontsize=15)
#plt.title(label)
#plt.savefig('plots/clf_cen_fit_'+str(M_host_min)+'-'+str(M_host_max)+'.png')
plt.show()
plt.close() 
"""


