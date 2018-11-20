#!/usr/bin/env python
# -*- coding: utf8 -*-

import os,sys
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

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


fig1 = plt.figure(figsize=(8,16))

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
plt.legend()
plt.grid()

# -------------------------------------------------------------------------------------
# Plot Ritzspectrum
# -------------------------------------------------------------------------------------
plt.subplot(212)
plt.scatter(lambdaDMD[1],[i/(2*np.pi) for i in lambdaDMD[0]], s=[40.*(i+.4) for i in amplog] , c=[40.*(i+.4) for i in amplog], marker='o')

plt.axhline(y=0.0, color='k', linestyle='--')


plt.xlabel('$\omega_i/2\pi$',fontsize=18)
plt.ylabel('$\omega_r$',fontsize=18)
plt.grid()
plt.show()

# # plt.xticks(fontsize=14)
# # plt.ylabel('Sound Pressure Level - SPL [dB]',fontsize=18)
# # plt.ylim(0,110)
# # plt.yticks(fontsize=14)

# plt.title('Cavity Sound Spectrum s',fontsize = 20)
# plt.legend()
# plt.grid()
# plt.show()
