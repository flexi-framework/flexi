#!/usr/bin/env python
# -*- coding: utf8 -*-

import os,sys
import argparse
import glob
import numpy as np
sys.ps1 = 'SOMETHING'
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# matplotlib.interactive(True)
parser = argparse.ArgumentParser(description='Plot DMD data 1.Eigenvalues, 2.Ritzspectrum')
parser.add_argument('-d', '--dmdFile', help='DMD File to be plotted')
args = parser.parse_args()

# -------------------------------------------------------------------------------------
# Data ReadIn
# -------------------------------------------------------------------------------------
dmdFile = open(args.dmdFile, "r") 
dmdData = dmdFile.readlines()
values = []
for line in dmdData[12:]:
  values.append(line.split())


alphaDMD  = [[float(i[0]) for i in values],[float(i[1]) for i in values]]
lambdaDMD = [[float(i[2]) for i in values],[float(i[3]) for i in values]]
sigmaDMD  = [[float(i[4]) for i in values],[float(i[5]) for i in values]]

amplog = [np.log10(np.sqrt(alphaDMD[0][i]*alphaDMD[0][i]+alphaDMD[1][i]*alphaDMD[1][i])) for i in range(len(alphaDMD[0]))]
maxA=max(amplog[2:])
minA=maxA-3               
for i in range(len(amplog)):
    amplog[i]=max(0.,min(1.,(amplog[i]-minA)/(maxA-minA)))



# fig1 = plt.figure(figsize=(32,16))

# -------------------------------------------------------------------------------------
# Plot Eigenvalues
# -------------------------------------------------------------------------------------
plt.subplot(211)
circle1=plt.Circle((0, 0), 1., color='k',fill=False)
plt.scatter(sigmaDMD[0],sigmaDMD[1], s=[40.*(i+.4) for i in amplog] , c=[40.*(i+.4) for i in amplog], marker='o')

plt.xlabel('$\sigma_r$',fontsize=18)
plt.ylabel('$\sigma_i$',fontsize=18)
plt.axis('equal')
plt.axis([-1.4, 1.4, -1.4, 1.4])
plt.gcf().gca().add_artist(circle1)
# plt.title('Eigenvalues ',fontsize = 20)
# plt.legend()
plt.grid()

# -------------------------------------------------------------------------------------
# Plot Ritzspectrum
# -------------------------------------------------------------------------------------
plt.subplot(212)
lambdaDMDreal = []
lambdaDMDimag = []
amplogPositiv =[]
for i in range(len(lambdaDMD[1])):
    if lambdaDMD[1][i] >= 0.00000:
        lambdaDMDreal.append(lambdaDMD[0][i])
        lambdaDMDimag.append(lambdaDMD[1][i])
        amplogPositiv.append(amplog[i])

plt.scatter([i/(2*np.pi) for i in lambdaDMDimag], lambdaDMDreal, s=[40.*(i+.4) for i in amplogPositiv] , c=[40.*(i+.4) for i in amplogPositiv], marker='o')
labels = ['Mode%d\n f=%.2f'%(i+1,lambdaDMDimag[i]/(2*np.pi)) for i in range(len(lambdaDMDreal))]
plt.axhline(y=0.0, color='k', linestyle='--')
# j=0
for label, x, y in zip(labels, [i/(2*np.pi) for i in lambdaDMDimag], lambdaDMDreal):
    # j+=1
    # i=(j)%2
    plt.annotate(
        label,
        xy=(x, y), xytext=(0, -20),
        textcoords='offset points', ha='center', va='top',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.1),
        arrowprops=dict(arrowstyle = '->', connectionstyle='arc3,rad=0'))

plt.xlabel('$\omega_i/2\pi$',fontsize=18)
plt.ylabel('$\omega_r$',fontsize=18)
# plt.xlim(0,6000)
# plt.ylim(-600,100)
plt.grid()
# plt.show()

def computeNRoomFreqs(n,c,geo):
    freq=[]
    for i in range(n):
        for j in range(0,2*n,2):
            for k in range(n):
                freqtmp=c/2.*np.sqrt((i/geo.x)**2.+((j+1)/(2*geo.y))**2.+(k/geo.z)**2.)
                freq.append([freqtmp,i,j+1,k])
    
    freq=sorted(freq, key=lambda x : x[0])
    
    return freq

# # def computeNRossiterModes(n,c,geo):

    # # return freq


# # def computeHelmholtzFreq(c,geo,neck):

    # # return freq



class room(object):
  def __init__(self,x,y,z):
      self.x = x
      self.y = y
      self.z = z

cavityGeo = room(0.025,0.05,0.03)

c=343.

freq = computeNRoomFreqs(4,c,cavityGeo)

for i in freq:
    plt.axvline(x=i[0], color='k', linestyle='--',label='Mode:'+str(i[1])+str(i[2])+str(i[3]))


# # plt.xticks(fontsize=14)
# # plt.ylabel('Sound Pressure Level - SPL [dB]',fontsize=18)
# # plt.ylim(0,110)
# # plt.yticks(fontsize=14)

# plt.title('Cavity Sound Spectrum s',fontsize = 20)
# plt.legend()
# plt.grid()
plt.show()
