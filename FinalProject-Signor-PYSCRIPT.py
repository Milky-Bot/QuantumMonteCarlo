#Short script which compiles and executes the fortran script
#"FinalProject-Signor-CODE.f90" which runs a Diffusion Monte Carlo Simulation
#to compute the ground state of many-body quantum systems.
#after the fortran program is completed, 1)the ground state energy is computed as
#the mean of the reference energy ("ref_energy.txt") and its error is the standard deviation
#2)the already computed ground state wavefunction is the spatial distribution of replicas ('wavefunction.txt')
#Comparsion via known results if available is delivered.
#Furthermore, the trace for the number of alive replicas is plotted, as diagnosis tool

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import subprocess
import os
import scipy.integrate as integrate
import pathlib
plt.rcParams.update({'font.size': 18})
sns.set_style('dark')


#%%Defining some exact ground states wavefunctions
def ha_gs(x):
    return np.pi**(-0.25)*np.exp(-0.5*x**2)
def morse_gs(x):
    return np.sqrt(2.)*np.exp(-np.exp(-x)-x/2.)
def Hydrogen_gs(x):
    return(2*np.exp(-x))


gs_wave_functions=[ha_gs,morse_gs,Hydrogen_gs]
#exact ground state energies, in order
#1:1D Harmonic Oscillator
#2: Morse potential
#3: Hidrogen atom
#4: H- ion
#5: H2+ ion
#6: H3+ ion
#7: H2 molecule
gs_energy=[0.5,-0.125,-0.5,-0.528,-0.597,-1.344,-1.161] 
#%%Running FORTRAN script
os.system("gfortran FinalProject-Signor-CODE.f90 -o qmc.out")
subprocess.call("./qmc.out")

#reading fortran script output
potential_choice=int(np.genfromtxt('potential_choice.txt')) # choice of potential
tau0=int(np.genfromtxt('time_steps.txt')) #number of time steps | fortran program input parameter
os.remove("potential_choice.txt");os.remove("time_steps.txt")

N=np.genfromtxt('N.txt') #number of alive replicas for timesteps >=tau0
Er=np.genfromtxt('ref_energy.txt') #reference energy
er = np.cumsum(Er)/np.arange(1,len(Er)+1)  #medium of reference energy


#%% Plots

#create a dedicated folder for plots
pathlib.Path('plots').mkdir(parents=True, exist_ok=True) 

#reference energy(time steps)
fig, ax = plt.subplots(1,figsize=(10,7))
ax.plot(np.arange(2*tau0),Er,c='red',label='Reference Energy')
ax.set_xlabel('\u03C4',fontsize=22)
ax.set_ylabel('$E_R$',fontsize=22)
ax.legend()
fig.tight_layout()
plt.savefig('plots/trace_ref_energy.png')

#number of alive replicas(time steps)
fig, ax = plt.subplots(1,figsize=(10,7))
ax.plot(np.arange(2*tau0),N,label='N')
ax.set_ylabel('N')

ax.set_xlabel('\u03C4')
fig.tight_layout()
plt.savefig('plots/alive_replicas.png')
print('Estimated E_0 = %0.3f'%np.mean(Er[int(len(Er)/2):]), 'Std = %0.3f'%np.std(Er[int(len(Er)/2):]) )

if potential_choice<=2:
    xmin=-10
    xmax=10
else:
    xmin=0
    xmax=10
if potential_choice<4:
    #mean reference energy(time_steps)
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(17,7))

    ax1.plot(np.arange(2*tau0),np.full(2*tau0,gs_energy[potential_choice-1]),label='$E_0^{ex}$')
    ax1.plot(np.arange(2*tau0),er,c='red',label='Mean reference energy')
    ax1.set_xlabel('\u03C4',fontsize=22)
    ax1.set_ylabel('$\overline{E_R}$',fontsize=22)
    ax1.legend()
   
    
    replicas_distribution=np.genfromtxt('wavefunction.txt') #gs wavefunction from simulation
    positions=replicas_distribution[:,0]
    wavefunction=replicas_distribution[:,1]

    x=np.linspace(xmin,xmax,1000)
    
    psi0=gs_wave_functions[potential_choice-1](x)
    psi0 /= integrate.quad(gs_wave_functions[potential_choice-1],xmin,xmax)[0]
    ax2.plot(x,psi0,label='Theoretical Result')
    ax2.set_xlim(xmin,xmax)    
    ax2.scatter(positions, wavefunction,s=4,c='red',label='Simulation')
    ax2.legend()
    ax2.set_xlabel('x',fontsize=22)
    ax2.set_ylabel('$\u03A6_0$',fontsize=22)
    fig.tight_layout()

    plt.savefig('plots/GS.png')


#for systems with 2 or more electrons
elif potential_choice >= 4:
    #mean reference energy(time_steps)
    fig, ax1= plt.subplots(1,1,figsize=(10,7))
    ax1.plot(np.arange(2*tau0),np.full(2*tau0,gs_energy[potential_choice-1]),label='$E_0^{ex}$')
    ax1.plot(np.arange(2*tau0),er,c='red',label='Mean reference energy')
    ax1.set_xlabel('\u03C4',fontsize=22)
    ax1.set_ylabel('$\overline{E_R}$',fontsize=22)
    ax1.legend()
    plt.savefig('plots/Gs-energy')
    if potential_choice>5:
        from mpl_toolkits.mplot3d import Axes3D
        spatial_coordinates=np.genfromtxt('positions.txt')
        r1=np.linalg.norm(spatial_coordinates[:,:3],axis=1) #first electron radius
        r2=np.linalg.norm(spatial_coordinates[:,3:],axis=1) #second electron radius

        
        fig = plt.figure(figsize=(17,7))
        #histoogram of 2 radii
        ax = fig.add_subplot(1, 2, 1)
        ax.hist2d(r1, r2, bins=200, range=((0,xmax),(0,xmax)), density=True)
        ax.set_xlabel('$R_1$')
        ax.set_ylabel('$R_2$')
        #3d distribution
        plt.rcParams.update({'font.size': 10})
        xs1 = spatial_coordinates[:,0]
        ys1 = spatial_coordinates[:,1]
        zs1 = spatial_coordinates[:,2]
        ax = fig.add_subplot(1,2,2,projection='3d')    
        ax.scatter(xs1, ys1, zs1,s=0.1,alpha=0.5,c='k')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        fig.tight_layout()
        plt.savefig('plots/spat_distr.png')
