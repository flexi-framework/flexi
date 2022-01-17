import argparse
import numpy as np
import h5py
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='ChannelFFT plotting routine')
optional = parser._action_groups.pop() 
required = parser.add_argument_group('required arguments')
required.add_argument('-p' ,'--projectName',help='Project name of the ChannelFFT output HDF5 files.',required=True)
required.add_argument('-t' ,'--time'       ,help='Time of the ChannelFFT output HDF5 files (arbitrary format float).',required=True)
optional.add_argument('-np','--noProfiles' ,help='Do no plot profiles of mean velocity and Reynolds stresses.', action='store_true')
optional.add_argument('-ns','--noSpectra'  ,help='Do not plot turbulent energy spectra.', action='store_true')
parser._action_groups.append(optional)

args = parser.parse_args()
time = str(format(float(args.time), '017.9f'))

if not args.noProfiles :
    FileName = args.projectName+"_MS_"+time+'.h5'
    DatasetName = 'MeanSquares'
    h5file = h5py.File(FileName, 'r')
    VarNames = h5file.attrs["VarNames"]
    meanSquares = np.transpose(np.array(h5file[DatasetName]))
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    for i in range(1,5):
        ax1.plot(meanSquares[0],meanSquares[i],label=VarNames[i])
    for i in range(7,8):
        ax2.plot(meanSquares[0],meanSquares[i], '--',label=VarNames[i])

    ax1.set_xlabel('$y^+$',fontsize = 16)
    ax1.set_ylabel('$\overline{u\'u\'}^+$,$\overline{v\'v\'}^+$,$\overline{w\'w\'}^+$,$\overline{u\'v\'}^+$',fontsize = 16)
    ax2.set_ylabel('$\overline{u}^+$',fontsize = 16)
    ax1.set_xlim(0,max(meanSquares[0]))
    ax1.legend(fontsize=12)
    ax2.legend(fontsize=12)
    plt.savefig(args.projectName+"_MS_"+time+'.png')
    plt.clf()

if not args.noSpectra:
    for dim in ["x","z"] :
        FileName = args.projectName+"_EnergySpectra_"+dim+"_"+time+'.h5'
        h5file = h5py.File(FileName, 'r')
        VarNames = h5file.attrs["VarNames"]
        for DatasetName in h5file.keys():
            energySpectra = np.transpose(np.array(h5file[DatasetName]))
            for i,Spec in enumerate(energySpectra[1:]):
                plt.loglog(energySpectra[0][1:],Spec[1:],label=VarNames[i+1])
            plt.xlabel('Wavenumber k',fontsize=16)
            plt.ylabel('E_uu_'+dim+',E_vv_'+dim+',E_ww_'+dim+',E_pp_'+dim,fontsize=16)
            plt.legend(fontsize=12)
            plt.title('yPlus='+DatasetName.split(' ')[-1],fontsize=16)
            plt.savefig(args.projectName+"_EnergySpectra_"+dim+"_yPlus"+DatasetName.split(' ')[-1]+"_"+time+".png")
            plt.clf()

