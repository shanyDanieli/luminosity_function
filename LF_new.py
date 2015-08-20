#!/usr/bin/env python                                                                                                                                                                 
###################                                                                                                                                                               
# Shany Danieli   #
# Yale University #                                                                                                                                                               
# 8/17/2015       #
###################


# Import python modules
from pylab import *
import numpy as np
from math import *
import matplotlib.pyplot as plt
from scipy import interpolate
import h5py

V = 250.0**3 # box volume

"""
importing data from files
"""
# only centrals
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
#data1_sat = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/splitted/sat_mag.txt",float)
#data2_sat = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/splitted/sat.txt",float)

# downsized sample
data1_sat = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/downsized/sat_mag_ds.txt",float)
data2_sat = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/downsized/sat_ds.txt",float)


abs_mag_sat = data1_sat[:,7]
Msun = 4.76
L_data_sat = ((Msun-abs_mag_sat)/2.5) # L is log(L/L_sun)
host_mass_sat = data2_sat[:,2]
upid_sat = data2_sat[:,1]


# Host mass bins
num_bins = 20
Mmin = 11.75
Mmax = 14.25
dlogM = (Mmax-Mmin)/num_bins
Mass_bins = np.arange(Mmin,Mmax+dlogM/2,dlogM)
Mass_centers = (Mass_bins+dlogM/2).tolist()
del Mass_centers[-1]

"""
#Mmin = 10.4
#Mmax = 15.0
#Mmin = 11.85
#Mmax = 14.15
dlogM = (Mmax-Mmin)/(num_bins-1)
Mass_centers = np.arange(Mmin,Mmax+dlogM/2,dlogM)
Mass_bins = (Mass_centers-dlogM/2).tolist()
Mass_bins.append(Mass_bins[-1]+dlogM)
"""

# Luminosity bins
L_min = 9.5 # 9.8	
L_max = 10.7	
dL = 0.15
L_bins = np.linspace(L_min,L_max,(L_max-L_min)/dL+1)
L_centers = [x+dL/2 for x in L_bins]
del L_centers[-1]




"""
calculates the satellites LF directly computed from the simulation 
"""

(counts, L_bins, patches) = hist(L_data_sat,L_bins) 
lumi_direct_sat = counts/(dL*V)
plt.close()



"""
calculates the CLF for satellite galaxies
"""
clf_sat = []
for x in range(num_bins):
	L_sat_binned = [L_data_sat[i] for i in range(len(L_data_sat)) if  Mass_bins[x]<log(host_mass_sat[i],10)<=Mass_bins[x+1]]
	id_sat_binned = [upid_sat[i] for i in range(len(upid_sat)) if  Mass_bins[x]<log(host_mass_sat[i],10)<=Mass_bins[x+1]]
	id_cen_binned = [upid_cen[i] for i in range(len(upid_cen)) if  Mass_bins[x]<log(host_mass_cen[i],10)<=Mass_bins[x+1]]
#	Nhosts = len(np.unique(id_sat_binned))
	Nhosts = len(np.unique(id_cen_binned))	
	(k, L_bins, patches) = hist(L_sat_binned,L_bins) 
	plt.close()
	clf_tmp = k/(dL*Nhosts)
	clf_sat.append(clf_tmp)
	


"""
calculates the halo mass function
"""
counts2 = np.histogram(np.log10(host_mass_cen), Mass_bins)[0]
plt.close()
dndlogM = counts2/(dlogM*V)

print "dndlogM"
print dndlogM
quit()

"""
calculated the LF from the CLF and HMF (satellites)
"""
lumi_from_clf_sat = sum([clf_sat[i]*dndlogM[i] for i in range(len(dndlogM))],axis=0)*dlogM

"""
plot both for a comparison
"""

plt.title('Satellite galaxies')
plt.scatter(L_centers,np.log10(lumi_from_clf_sat),color='b' ,label='from clf')
plt.scatter(L_centers,np.log10(lumi_direct_sat),color='r' ,label='direct')
plt.legend()
plt.show()















"""
Everything for central galaxies
"""

"""
calculates the satellites LF directly computed from the simulation 
"""
(counts, L_bins, patches) = hist(L_data_cen,L_bins) 
lumi_direct_cen = counts/(dL*V)
plt.close()


"""
calculates the CLF for central galaxies
"""
clf_cen = []
for x in range(num_bins):
	L_cen_binned = [L_data_cen[i] for i in range(len(L_data_cen)) if  Mass_bins[x]<log(host_mass_cen[i],10)<=Mass_bins[x+1]]
	id_cen_binned = [upid_cen[i] for i in range(len(upid_cen)) if  Mass_bins[x]<log(host_mass_cen[i],10)<=Mass_bins[x+1]]
	Nhosts = len(np.unique(id_cen_binned))
	(k, L_bins, patches) = hist(L_cen_binned,L_bins) 
	plt.close()
	clf_tmp = k/(dL*Nhosts)
	clf_cen.append(clf_tmp)
	

"""
calculated the LF from the CLF and HMF (satellites)
"""
lumi_from_clf_cen = sum([clf_cen[i]*dndlogM[i] for i in range(len(dndlogM))],axis=0)*dlogM


"""
plot both for a comparison
"""
plt.title('Central galaxies')
plt.scatter(L_centers,np.log10(lumi_from_clf_cen),color='b' ,label='from clf')
plt.scatter(L_centers,np.log10(lumi_direct_cen),color='r' ,label='direct')
plt.legend()
plt.show()


print "CLF:"
print "satellites:"
print clf_sat
print "cenetrals:"
print clf_cen



print "Computed directly"
print "lumi_direct_sat:"
print lumi_direct_sat
print "lumi_direct_cen:"
print lumi_direct_cen



print "Computed by integrating the CLF and HMF"
print "lumi_from_clf_sat:"
print lumi_from_clf_sat
print "lumi_from_clf_cen:"
print lumi_from_clf_cen


print "diff:"
print 100*(lumi_direct_cen-lumi_from_clf_cen)/lumi_direct_cen


print "\n\n"
print "************************* DONE *************************"







