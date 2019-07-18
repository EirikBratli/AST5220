import numpy as np
import matplotlib.pyplot as plt


def Plotting(file, k_file, x_file, name):
    data = np.loadtxt(file)
    ks, ckH0, k_i = np.loadtxt(k_file, unpack=True)
    x, ka, kb = np.loadtxt(x_file, unpack=True)
    t2 = np.loadtxt('theta2.dat', unpack=True)
    color = ['b', 'g', 'r', 'c', 'm', 'y']
    
    print(np.shape(data))
   
    if name == 'v' or name == 'v_b':
        ylab = name
    else:
        
        ylab = ('\%s'%name)
    print(ylab)  
  
    for i in range(0,6):
        print(i, '%g'%ckH0[i])
        
        plt.figure('log %s'%name)
        plt.semilogy(x[:], data[:,i], c=color[i], label=r'$ck/H_0$=%g'%ckH0[i])
        plt.xlabel(r'$x$', size=15)
        plt.ylabel(r'$%s$'%ylab, size=15)
        plt.legend(loc='best', fontsize=12)
        plt.grid(True)
        plt.savefig('Plots_mk3/logplot_%s.png'%name)

        plt.figure('num $%s$'%name)
        plt.plot(x[:], data[:,i], c=color[i], label=r'$ck/H_0$=%g'%ckH0[i])
        plt.xlabel(r'$x$', size=15)
        plt.ylabel(r'$%s$'%ylab, size=15)
        plt.legend(loc='best', fontsize=12)
        plt.grid(True)
        plt.savefig('Plots_mk3/plot_%s.png'%name)
    
    plt.show()
    

Plotting('data_mk3/phi.dat', 'data_mk3/k_values.dat', 'data_mk3/x_t.dat', 'Phi')
Plotting('data_mk3/psi.dat', 'data_mk3/k_values.dat', 'data_mk3/x_t.dat', 'Psi')
Plotting('data_mk3/delta.dat', 'data_mk3/k_values.dat', 'data_mk3/x_t.dat', 'delta')
Plotting('data_mk3/delta_b.dat', 'data_mk3/k_values.dat', 'data_mk3/x_t.dat', 'delta_b')
Plotting('data_mk3/v.dat', 'data_mk3/k_values.dat', 'data_mk3/x_t.dat', 'v')
Plotting('data_mk3/v_b.dat', 'data_mk3/k_values.dat', 'data_mk3/x_t.dat', 'v_b')
Plotting('data_mk3/Theta0.dat', 'data_mk3/k_values.dat', 'data_mk3/x_t.dat', 'Theta_0')
Plotting('data_mk3/Theta1.dat', 'data_mk3/k_values.dat', 'data_mk3/x_t.dat', 'Theta_1')

