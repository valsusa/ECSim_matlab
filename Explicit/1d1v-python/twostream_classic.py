""" 
    1D electrostatic particle-in-cell solver for electron two-stream instability.

    Les Houches by G. Lapenta.

    V. Olshevsky: sya@mao.kiev.ua

"""
import os, time
start_time = time.clock()
import numpy as np
import pylab as plt
import matplotlib.patches as mpatches
from scipy import sparse
from scipy.sparse import linalg

# Output folder
path = 'C:/SyA/students/Leuven/python/output'
if not os.path.exists(path):
    os.makedirs(path)

# Set plotting parameters
params = {'axes.labelsize': 'large',
              'xtick.labelsize': 'medium',
              'ytick.labelsize': 'medium',
              'font.size': 15,
              'font.family': 'sans-serif',
              'text.usetex': False,
              'mathtext.fontset': 'stixsans',}
plt.rcParams.update(params)
## Switch on interactive plotting mode
plt.ion()

# Simulation parameters
L = 20*np.pi #20*np.pi # Domain size
DT = 0.005 # Time step
NT = 50000  # Number of time steps
TOut = round(NT/500) # Output period
verbose = True
NG = 320  # Number of grid cells
N = NG * 20 # Number of particles
WP = 1. # Plasma frequency
QM = -1. # Charge/mass ratio
V0 = 0.9 # Stream velocity
VT = 0.0000001 # Thermal speed

# perturbation 
XP1 = 1.0 
mode = 1

Q = WP**2 / (QM*N/L)  # rho0*L/N: charge carried by a single particle?
rho_back = -Q*N/L  # Background charge density?
dx = L / NG # Grid step

# Auxilliary vectors
p = np.concatenate([np.arange(N), np.arange(N)])  # Some indices up to N
Poisson = sparse.spdiags(([1, -2, 1] * np.ones((1, NG-1), dtype=int).T).T, \
                         [-1, 0, 1], NG-1, NG-1)
Poisson = Poisson.tocsc()

# Cell center coordinates
xg = np.linspace(0, L-dx, NG) + dx/2

# electrons
xp = np.linspace(0, L-L/N, N).T   # Particle positions
vp = VT * (1 - VT**2)**(-0.5) * np.random.randn(N) # Particle momentum, initially Maxwellian
pm = np.arange(N)
pm = 1 - 2 * np.mod(pm+1, 2)
vp += pm * (V0 * (1 - V0**2)**(-0.5)) # Momentum + stream velocity

# Add electron perturbation to excite the desired mode
xp += XP1 * (L/N) * np.sin(2 * np.pi * xp / L * mode)
xp[np.where(xp < 0)] += L
xp[np.where(xp >= L)] -= L

histEnergy, histPotE, histKinE, histMomentum, t = [], [], [], [], []

if verbose:
    plt.figure(1, figsize=(16,9))

# Main cycle
for it in xrange(NT+1):      
    # update particle position xp
    xp += vp * DT
    # Periodic boundary condition
    xp[np.where(xp < 0)] += L
    xp[np.where(xp >= L)] -= L

    # Project particles->grid
    g1 = np.floor(xp/dx - 0.5)
    g = np.concatenate((g1, g1+1))
    fraz1 = 1 - np.abs(xp/dx - g1 - 0.5)
    fraz = np.concatenate((fraz1, 1-fraz1))
    g[np.where(g < 0)] += NG
    g[np.where(g > NG-1)] -= NG

    mat = sparse.csc_matrix((fraz, (p, g)), shape=(N, NG))
    rho = Q / dx * mat.toarray().sum(axis=0) + rho_back 

    # Compute electric field potential
    Phi = linalg.spsolve(Poisson, -dx**2 * rho[0:NG-1])
    Phi = np.concatenate((Phi,[0]))
    # Electric field on the grid
    Eg = (np.roll(Phi, 1) - np.roll(Phi, -1)) / (2*dx)
   
    # interpolation grid->particle and velocity update
    vp += mat * QM * Eg * DT
   
    Etot = 0.5 * (Eg**2).sum() * dx
    histEnergy.append(Etot)
    histPotE.append(0.5 * (Eg**2).sum() * dx)
    histKinE.append(0.5 * Q/QM * (vp**2).sum())
    histMomentum.append(Q/QM * vp.sum())
    t.append(it*DT)
   
    if (np.mod(it, TOut) == 0) and verbose:
        # Phase space        
        plt.clf()
        plt.subplot(2, 2, 1)
        plt.scatter(xp[0:-1:2], vp[0:-1:2], s=0.5, marker='.', color='blue')
        plt.scatter(xp[1:-1:2], vp[1:-1:2], s=0.5, marker='.', color='red')
        plt.xlim(0, L)
        plt.ylim(-7, 7)
        plt.xlabel('X')
        plt.ylabel('P')
        plt.legend((mpatches.Patch(color='w'), ), (r'$\omega_{pe}t=$' + str(DT*it), ), loc=1, frameon=False)
        
        # Electric field
        plt.subplot(2, 2, 2)
        plt.xlim(0, L)
        plt.xlabel('X')
        plt.plot(xg, Eg, label='E', linewidth=2)
        plt.legend(loc=1)
        
        # Energies
        plt.subplot(2, 2, 3)
        plt.xlim(0, NT*DT)
        plt.xlabel('time')
        plt.yscale('log')
        plt.plot(t, histEnergy, label='Total Energy', linewidth=2)
        plt.plot(t, histPotE, label='Potential', linewidth=2)
        plt.plot(t, histKinE, label='Kinetic', linewidth=2)
        plt.legend(loc=4)
        
        # Momentum
        plt.subplot(2, 2, 4)
        plt.xlim(0, NT*DT)
        plt.xlabel('time')
        plt.plot(t, histMomentum, label='Momentum', linewidth=2)
        plt.legend(loc=1)
        plt.pause(0.01)
        print it
        #plt.savefig(os.path.join(path, 'twostream%3.3i' % (it/TOut,) + '.png'))
        

print 'Time elapsed: ', time.clock() - start_time