import matplotlib.pyplot as plt
import numpy
import pandas

advvel   = 1.   # advection velocity
t        = 5.   # time stamp of visualization

# plot numerical and analytical solution
def PlotSolution(filename,ax,plotlabel,fn_A):
    df = pandas.read_csv(filename)
    ax.plot(df['Points_0'], fn_A(df['Points_0']),label='exact')
    ax.plot(df['Points_0'], df['Solution']      ,label=plotlabel)
    ax.legend(loc='upper right')
    ax.set_xlabel("x")
    ax.set_ylabel("A")

# well-resolved simulation
filename  = 'LinAdv_t5_WellResolved_N4.csv'
plotlabel = 'Solution N=4'
omega     = 1.96349540849
fn_A      = lambda x: numpy.cos(omega/advvel*(x-advvel*t))   # define analytical solution for amplitude (real part of solution)
fig,axs   = plt.subplots(1,1)
fig.set_size_inches(8,4.5)
PlotSolution(filename,axs,plotlabel,fn_A)
plt.savefig(filename[:-4]+'.jpg', bbox_inches='tight')

# under-resolved simulations
filename  = ['LinAdv_t5_UnderResolved_N2.csv','LinAdv_t5_UnderResolved_N4.csv','LinAdv_t5_UnderResolved_N6.csv','LinAdv_t5_UnderResolved_N11.csv']
plotlabel = ['Solution N=2','Solution N=4','Solution N=6','Solution N=11']
omega     = [2.4,4.0,5.6,9.6]
fig,axs   = plt.subplots(2,2)
fig.set_size_inches(16,9)
for n in range(len(filename)):
    fn_A = lambda x: numpy.cos(omega[n]/advvel*(x-advvel*t))   # redefine analytical solution for amplitude (real part of solution)
    PlotSolution(filename[n],numpy.ravel(axs)[n],plotlabel[n],fn_A)
plt.savefig('LinAdv_t5_UnderResolved.jpg', bbox_inches='tight')