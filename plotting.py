import numpy as np
import matplotlib.pyplot as plt

def cmb_data():
    """
    Imported file contains: x_eta, z, eta, eta2, H(x_eta), H'(x_eta), dH'(x_eta), 
    
    """

    data = np.loadtxt('cmb_data.dat')
    names = [r'$x_{\eta}$',r'$z$',r'$\eta$',r'$H(x_{\eta})$',r'$\mathcal{H}(x_{\eta})$', r'd$\mathcal{H}(x_{\eta})$']
    
    namez = [r'$x_{\eta} = \ln{a}$',r'$z$',r'$\eta$',r'$H(z)$',r'$\mathcal{H}(z)$',r'd$\mathcal{H}(z)$']
    save = ['x', 'z', 'eta1', 'H', 'H_prime', 'dH_prime']
    print(np.shape(data))

    plt.figure('H(z)')
    plt.loglog(data[:,1], data[:,3], '-r')
    plt.xlabel(r'$z$', size=15)
    plt.ylabel(r'$H(z)$', size=15)
    plt.grid(True)
    plt.gca().invert_xaxis()
    plt.savefig('Plots_mk1/H_z.pdf')

    plt.figure('H(x)')
    plt.semilogy(data[:,0], data[:,3], '-b')
    plt.xlabel(r'$x$', size=15)
    plt.ylabel(r'$H(x)$', size=15)
    plt.grid(True)
    plt.savefig('Plots_mk1/H_x.pdf')

    plt.figure('eta(x)')
    plt.semilogy(data[:,0], data[:,2], '-b')
    plt.xlabel(r'$x$', size=15)
    plt.ylabel(r'$\eta(x)$', size=15)
    plt.grid(True)
    plt.savefig('Plots_mk1/eta_x.pdf')


def eta_xt():
    """
    Imported file contains: x_t, z, eta(x_t)
    """
    
    data = np.loadtxt('eta.dat')
    print(np.shape(data[:,0]))
    
    plt.figure('eta(x)')
    plt.semilogy(data[:,0], data[:,2], '-g')
    plt.xlabel(r'$x$', size=15)
    plt.ylabel(r'$\eta(x)$', size=15)
    plt.grid(True)
    #plt.savefig('Plots_mk1/eta_x.pdf')


def densities():
    """
    Plotting the densities as function of red shift or x
    """
    
    data = np.loadtxt('densities.dat')
    Omega_i = ['-','-', r'$\Omega_r$', r'$\Omega_m$', r'$\Omega_b$', r'$\Omega_{\Lambda}$']
    color = ['m','m', 'b', 'r', 'g', 'k']
    print(np.shape(data[:,:]))
    plt.figure('densities')
    for i in range(2, len(data[0,:])):
        plt.semilogx(data[:,0], data[:,i], c=color[i], label='%s'%Omega_i[i])
        
    plt.xlabel(r'$z$', size=15)
    plt.ylabel(r'$\Omega(z)$ / $\sum{\Omega_i}$', size=15)
    plt.legend(loc='best', fontsize=12)
    plt.gca().invert_xaxis()
    plt.grid(True)
    plt.savefig('Plots_mk1/densities.pdf')

    plt.figure('densities-x')
    for i in range(2, len(data[0,:])):
        plt.plot(data[:,1], data[:,i], c=color[i], label='%s'%Omega_i[i])
        
    plt.xlabel(r'$x$', size=15)
    plt.ylabel(r'$\Omega(x)$ / $\sum{\Omega_i}$', size=15)
    plt.legend(loc='best', fontsize=12)
    #plt.gca().invert_xaxis()
    plt.grid(True)
    plt.savefig('Plots_mk1/densities_x.pdf')


"""
Function calls
"""

#cmb_data()
#eta_xt()
densities()
plt.show()
