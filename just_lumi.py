#!/usr/bin/env python                                                                                                                                                                 
###################                                                                                                                                                               
# Shany Danieli   #
# Yale University #                                                                                                                                                               
# 7/27/2014       #
###################

# Import python modules
from pylab import *
from numpy import *
from math import *
import matplotlib.pyplot as plt


# Read in the data from the file "conditional_abundance_matching_Mr_gr_model.dat"
data = loadtxt("../../Data/erased_assembly_bias_Mr_gr_model.dat",float) #abundance matching mock with erased assembly bias
#data = loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/erased_assembly_bias_Mr_gr_model_downsized.dat",float) #abundance matching mock with erased assembly bias, downsized to Mhost=[11.75-14.25]
#data = loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/erased_assembly_bias_Mr_gr_model_downsized2.dat",float) #abundance matching mock with erased assembly bias, downsized to Mhost=[11.75-14.25]

data2 = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock.dat",float) # no assembly bias
#data2 = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock_downsized.dat",float) # no assembly bias, downsized to Mhost=[11.75-14.25]
#data2 = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock_downsized2.dat",float) # no assembly bias, downsized to Mhost=[11.75-14.25]
halo_ID = data2[:,0] # Halo ID
UPID = data2[:,1] # if UPID == -1 -> central, if UPID > 0 -> satellite
Mhost = data2[:,2] # Host mass in units of (Msun/h)

print len(data)
print len(data2)

absMag = data[:,7]
pos = data[:,range(1,4)] # the position of each galaxy in the box
N = len(absMag)

Vbox = 15625000 # The volume of the box
delta_M = 0.15 # size of the bins


# Calculation of the luminosity function
bins = [x*0.15-22 for x in range(0,21)]
(n, bins, patches) = hist(absMag,bins,color='w') # 0.15 bins from -19 to -22
Nbins = len(bins)
magnitude = [round(x-delta_M/2,3) for x in bins]
del magnitude[0]
plt.close()

#dndMag = n.tolist()
dndMag = [x/(Vbox*delta_M) for x in n] # the luminosity function



data = loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/erased_assembly_bias_Mr_gr_model_downsized3.dat",float) #abundance matching mock with erased assembly bias, downsized to Mhost=[11.75-14.25]
data2 = np.loadtxt("/Users/si224/Documents/2014/vanDenBosch/Data/noab_mock_downsized3.dat",float) # no assembly bias, downsized to Mhost=[11.75-14.25]
halo_ID = data2[:,0] # Halo ID
UPID = data2[:,1] # if UPID == -1 -> central, if UPID > 0 -> satellite
Mhost = data2[:,2] # Host mass in units of (Msun/h)
absMag = data[:,7]
pos = data[:,range(1,4)] # the position of each galaxy in the box
N = len(absMag)
bins = [x*0.15-22 for x in range(0,21)]
(n2, bins, patches) = hist(absMag,bins,color='w') # 0.15 bins from -19 to -22
magnitude = [round(x-delta_M/2,3) for x in bins]
del magnitude[0]
plt.close()

dndMag2 = [x/(Vbox*delta_M) for x in n2] # the luminosity function


plt.semilogy(magnitude,dndMag,marker='o',linestyle='none',c='r',label='full')
plt.semilogy(magnitude,dndMag2,marker='o',linestyle='none',c='b',label='downsized')
plt.legend()
plt.xlabel(r'$M_r-5\logh$',fontsize=15)
plt.ylabel(r'$\log(\Phi[h^{3}Mpc^{-3}mag^{-1}])$',fontsize=15)
plt.gca().invert_xaxis()
plt.show()





"""
# plot the luminosity function
# with the analytical calculation from the CLF and hmf
errors_log = [x/y for x,y in zip(errors,dndMag)]
plt.figure()
plt.subplot(2,1,1)
plt.errorbar(magnitude,lumi_log, yerr=errors_log, capsize=4, ls='none',color = 'k', elinewidth=1,marker='x',label='Bolshoi')
plt.scatter(magnitude,[log(x,10) for x in lumi_model],color = 'b',marker='o',label='Model(Tinker)')
plt.title(r'$\Phi(L)$')
plt.xlabel(r'$M_r-5\logh$',fontsize=15)
plt.ylabel(r'$\log(\Phi[h^{3}Mpc^{-3}mag^{-1}])$',fontsize=15)
plt.gca().invert_xaxis()
plt.legend()

plt.subplot(2,1,2)
plt.plot(magnitude,diff,linestyle='None',marker='o',markersize=5,markeredgecolor='k', markerfacecolor='m')
plt.xlabel(r'$M_r-5\logh$',fontsize=15)
plt.ylabel(r'$(\Phi_{Bolshoi}-\Phi_{Model})/\Phi_{Bolshoi})$',fontsize=15)
plt.ylim(-0.4,0.4)
plt.axhline(linewidth=0.5, color='k')
plt.axhspan(-0.05, 0.05, facecolor='0.5',linewidth=0)
plt.axhspan(-0.05, -0.1, facecolor='0.7',linewidth=0)
plt.axhspan(0.05, 0.1, facecolor='0.7',linewidth=0)
plt.axhspan(0.1, 0.2, facecolor='0.9',linewidth=0)
plt.axhspan(-0.1, -0.2, facecolor='0.9',linewidth=0)
plt.gca().invert_xaxis()
plt.savefig('directVsCLF_tinker_'+str(num_bins)+'.png')
"""

"""
plt.subplot(2,1,2)
plt.errorbar(magnitude,lumi_log, yerr=errors_log, capsize=4, ls='none',color = 'k', elinewidth=1,marker='x')
plt.scatter(magnitude,[log(x,10) for x in lumi_data],color = 'b',marker='o',label='Mock')
plt.xlabel(r'$M_r-5\logh$',fontsize=15)
plt.ylabel(r'$\log(\Phi[h^{3}Mpc^{-3}mag^{-1}])$',fontsize=15)
plt.gca().invert_xaxis()
plt.legend()
"""

#plt.savefig('directVsCLF'+str(num_bins)+'.png')
plt.show()
plt.close()


