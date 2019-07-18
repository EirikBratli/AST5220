import numpy as np
import matplotlib.pyplot as plt


z, x , X_e, n_e, tau, tau2, tau22, g, g2, g22 = np.loadtxt('electron_density.dat', unpack=True)

def plot_Xe():
    """
    Function for plotting X_e as funciton of redshift
    """

    print('Plotting fractional electron density')

    plt.figure('Electron density')
    plt.semilogy(z, X_e, '-b')
    plt.xlabel(r'$z$', size=14)
    plt.ylabel(r'$X_e$', size=14)
    plt.grid(True)
    plt.xlim(0, 1800)
    plt.gca().invert_xaxis()
    #plt.savefig('Plots_mk2/electron_density.png')
    plt.show()


def plot_tau():
    print('Plotting the optical depth')

    plt.figure('Optical depth')
    plt.semilogy(x[:-10], tau[:-10], '-b', label=r'$\tau$')
    plt.semilogy(x[:-10], np.abs(tau2[:-10]), '-r', label=r'|d$\tau$|')
    plt.semilogy(x[:-10], np.abs(tau22[:-10]), '-g', label=r'|d$^2 \tau$|')
    plt.xlabel(r'$x$', size=14)
    plt.ylabel(r'Optical depth, $\tau$, |d$\tau$|', size=14)
    plt.grid(True)
    plt.legend(loc=1, fontsize=12)
    #plt.savefig('Plots_mk2/optical_depth.png')
    plt.show()

def plot_g():

    print('Plotting visibilty function')

    i_max = np.where(g==np.max(g))
    print g[i_max], x[i_max], np.exp(-x[i_max])-1
    print tau[i_max], X_e[i_max]

    plt.figure('Visibility function')
    plt.plot(x, g, '-b', label=r'$\tilde{g}$')
    plt.plot(x, g2/10, '-r', label=r'd$\tilde{g}/10$')
    plt.plot(x, g22/300, '-g', label=r'$d^2 \tilde{g}/300$')
    plt.xlabel(r'$x$', size=14)
    plt.ylabel(r'Visibility function, $\tilde{g}$', size=14)
    plt.legend(loc=1, fontsize=12)
    plt.xlim(-8,-6)
    plt.grid(True)
    plt.savefig('visibility_function2.png')
    plt.show()


#plot_Xe()
#plot_tau()
plot_g()
