#!/usr/bin/python


# Import from modules numpy and pylab
from pylab import *
from numpy import *
from matplotlib import pyplot 

left = [-22,-21,-20,-19.5]
lumiPerVol = [0.001070464,0.005077056,0.0043024,0.005184576]
lumiPerVolPerBin = [0.001070464,0.005077056,0.0086048,0.010369152]
yError = [0.00024076825518327787,0.0011541122555314973,0.00097516825411823156,0.0011697282556217918]
yErrorPerBin = [0.00024076825518327787,0.0011541122555314973, 0.00195033650823646312, 0.0023394565112435836]
width = [1,1,0.5,0.5]

#bar(left, lumiPerVol,width, yerr = yError, color='w',linewidth=1.5,error_kw=dict(elinewidth=4,ecolor='red'))
bar(left, lumiPerVolPerBin,width, yerr = yError, color='w',linewidth=1.5,error_kw=dict(elinewidth=4,ecolor='red'))
xlabel(r'$M_r-5logh$',fontsize=15)
ylabel(r'$N_{gal}[h^3Mpc^-3]$',fontsize=15)
xlim(-18,-23)
title('The luminosity function per unit volume',fontsize=15)
show()
