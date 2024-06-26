'''Example demonstrating Sedov solvers. Reproduces plots from Kamm & Timmes,
"On Efficient Generation of Numerically Robust Sedov Solutions," LA-UR-07-2849

Uses Doebling and (if available) Timmes Sedov solvers from ExactPack 2023 release at github.com.

ExactPack will need to be installed locally for this code to run.
'''


# import standard Python packages
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib import rc
rc('font', size=14)

# import ExactPack solvers
from exactpack.solvers.sedov.doebling import Sedov as SedovDoebling

# pyplot default settings
# rc('font', **{'family': 'serif', 'serif': ['Computer Modern'], 'size': 16})
# rc('grid', c='0.5', ls='-', lw=0.5)

# set domain variables for plots
npts = 1001
xl = 0.
xr = 1.2
rvec = np.linspace(0.0, 1.2, npts)
t = 1.0

#
# Figure 8doebling: Standard test cases, Doebling Solver
#

#solver_doebling_pla = SedovDoebling(geometry=1, eblast=0.0673185,
#                                    gamma=1.4, omega=0.)
#solution_doebling_pla = solver_doebling_pla(r=rvec, t=t)

solver_doebling_cyl = SedovDoebling(geometry=2, eblast=0.311357,
                                    gamma=1.4, omega=0.)
solution_doebling_cyl = solver_doebling_cyl(r=rvec, t=t)

rho = solution_doebling_cyl.density
V = solution_doebling_cyl.velocity
eint = solution_doebling_cyl.specific_internal_energy
P = solution_doebling_cyl.pressure
mom = solution_doebling_cyl.density*V
Cspd = solution_doebling_cyl.sound_speed
Erho = P/(0.4) + rho*(V*V)/2.
Ma = V/Cspd

# Write to output
filename = 'none'
icase = 'sedov'
vartype = 'CV'
gamma = 7./5.

if(filename == 'none'):
    filename = 'exact-{:s}-{:s}.dat'.format(vartype, icase)
    sys.stdout.write('  Writing to file: {:s}\n'.format(filename))
    with open(filename,'w+') as file:
        file.write('# Exact Riemann solution for problem: {:s}\n'.format(icase))
        file.write('# Discretization: {:s}\n'.format(str(npts)))
        file.write('# Domain bounds: [{:<.5e},{:.5e}]\n'.format(xl,xr))
        file.write('# Gamma: {:.5e}\n'.format(gamma))
        file.write('# Time: {:.5e}\n'.format(t))
        #file.write('# Iterface position: {:.5e}\n'.format(x0))
        #file.write('# Left  state [dens,velx,pres]: [{:<.5e}, {:.5e}, {:.5e}]\n'.format(stateL[0],stateL[1],stateL[2]))
        #file.write('# Right state [dens,velx,pres]: [{:<.5e}, {:.5e}, {:.5e}]\n'.format(stateR[0],stateR[1],stateR[2]))
        if(vartype == 'CV'):
            file.write('{:>20} {:>20} {:>20} {:>20} {:>20}\n'.format('x','rho','|m_i|','Erho','Ma'))
            for i in range(0,npts):                                                                                                                                       
                file.write('{:20.12e} {:20.12e} {:20.12e} {:20.12e} {:20.12e}\n'.format(rvec[i],rho[i],mom[i],Erho[i],Cspd[i]))
        else:
            file.write('{:>20} {:>20} {:>20} {:>20} {:>20} {:>20}\n'.format('x','rho','P','|v_i|','eint','Cspd'))
            for i in range(0,npts):
                file.write('{:20.12e} {:20.12e} {:20.12e} {:20.12e} {:20.12e} {:20.12e}\n'.format(rvec[i],rho[i],P[i],V[i],eint[i],Cspd[i]))

filename = 'none'
vartype = 'PV'
if(filename == 'none'):
    filename = 'exact-{:s}-{:s}.dat'.format(vartype, icase)
    sys.stdout.write('  Writing to file: {:s}\n'.format(filename))
    with open(filename,'w+') as file:
        file.write('# Exact Riemann solution for problem: {:s}\n'.format(icase))
        file.write('# Discretization: {:s}\n'.format(str(npts)))
        file.write('# Domain bounds: [{:<.5e},{:.5e}]\n'.format(xl,xr))
        file.write('# Gamma: {:.5e}\n'.format(gamma))
        file.write('# Time: {:.5e}\n'.format(t))
        #file.write('# Iterface position: {:.5e}\n'.format(x0))
        #file.write('# Left  state [dens,velx,pres]: [{:<.5e}, {:.5e}, {:.5e}]\n'.format(stateL[0],stateL[1],stateL[2]))
        #file.write('# Right state [dens,velx,pres]: [{:<.5e}, {:.5e}, {:.5e}]\n'.format(stateR[0],stateR[1],stateR[2]))
        if(vartype == 'CV'):
            file.write('{:>20} {:>20} {:>20} {:>20} {:>20}\n'.format('x','rho','|m_i|','Erho','Cspd'))
            for i in range(0,npts):                                                                                                                                       
                file.write('{:20.12e} {:20.12e} {:20.12e} {:20.12e} {:20.12e}\n'.format(rvec[i],rho[i],mom[i],Erho[i],Cspd[i]))
        else:
            file.write('{:>20} {:>20} {:>20} {:>20} {:>20} {:>20}\n'.format('x','rho','P','|v_i|','eint','Cspd'))
            for i in range(0,npts):
                file.write('{:20.12e} {:20.12e} {:20.12e} {:20.12e} {:20.12e} {:20.12e}\n'.format(rvec[i],rho[i],P[i],V[i],eint[i],Cspd[i]))

#print("solution: ", Cspd)
#solver_doebling_sph = SedovDoebling(geometry=3, eblast=0.851072,
#                                    gamma=1.4, omega=0.)
#solution_doebling_sph = solver_doebling_sph(r=rvec, t=t)

fig = plt.figure(figsize=(10, 7))
plt.suptitle('''Sedov solutions for $\gamma=1.4$, standard cases, Doebling solver.
    Compare to Fig. 8 from Kamm & Timmes 2007''')

plt.subplot(221)
#solution_doebling_pla.plot('density')
solution_doebling_cyl.plot('density')
#solution_doebling_sph.plot('density')
plt.xlim(0.0, 1.2)
plt.ylim(0.0, 6.5)
plt.xlabel('Position (cm)')
plt.ylabel('Density (g/cc)')
plt.grid(True)
L = plt.legend(loc='upper left', bbox_to_anchor=(0.25, 1.4), ncol=3,
               fancybox=True, shadow=True)
#L.get_texts()[0].set_text('planar')
#L.get_texts()[1].set_text('cylindrical')
L.get_texts()[0].set_text('cylindrical')
#L.get_texts()[2].set_text('spherical')

plt.subplot(222)
#solution_doebling_pla.plot('velocity')
solution_doebling_cyl.plot('velocity')
#solution_doebling_sph.plot('velocity')
plt.xlim(0.0, 1.2)
plt.ylim(0.0, 0.4)
plt.xlabel('Position (cm)')
plt.ylabel('Particle velocity (cm/s)')
plt.grid(True)
plt.gca().legend().set_visible(False)

plt.subplot(223)
#solution_doebling_pla.plot('specific_internal_energy')
solution_doebling_cyl.plot('specific_internal_energy')
#solution_doebling_sph.plot('specific_internal_energy')
plt.xlim(0.0, 1.2)
plt.ylim(1.e-2, 1.e5)
plt.xlabel('Position (cm)')
plt.ylabel('Specific internal energy (erg/g)')
plt.grid(True)
plt.gca().set_yscale('log', nonpositive='clip')
plt.gca().legend().set_visible(False)

plt.subplot(224)
#solution_doebling_pla.plot('pressure')
solution_doebling_cyl.plot('pressure')
#solution_doebling_sph.plot('pressure')
plt.xlim(0.0, 1.2)
plt.ylim(0.0, 0.15)
plt.xlabel('Position (cm)')
plt.ylabel('Pressure (erg/cc)')
plt.grid(True)
plt.gca().legend().set_visible(False)

plt.tight_layout()
fig.subplots_adjust(top=0.8)  # Makes room for suptitle
#plt.savefig('fig08doebling.pdf')
plt.show()

plt.close()
