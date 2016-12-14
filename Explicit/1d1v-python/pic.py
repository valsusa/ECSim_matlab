"""1D particle-in-cell solver for two-stream instability.

A translation of pic.m MATLAB routine by G. Lapenta.

V. Olshevsky: sya@mao.kiev.ua

"""
import time
start_time = time.clock()
import numpy as np
import pylab as plt
from scipy import sparse
from scipy.sparse import linalg

# Set plotting parameters
params = {'axes.labelsize': 'large',
              'xtick.labelsize': 'medium',
              'ytick.labelsize': 'medium',
              'font.size': 10,
              'font.family': 'sans-serif',
              'text.usetex': False,
              'mathtext.fontset': 'stixsans',}
plt.rcParams.update(params)
## Switch on interactive plotting mode
plt.ion()


L = 4*np.pi # Domain size
DT = 0.005 # Time step
NT = 20000  # Number of time steps
TOut = round(NT/100) # Output period
verbose = True
NG = 80 #320  # Number of grid cells
N = NG * 100 #100 # Number of particles
WP = 1.
QM = -1. # Charge/mass?
V0 = 0.5 # Stream velocity?
VT = 0.001 # Thermal speed?
# perturbation 
XP1 = 1.0 
V1 = 0.0
mode = 2

Q = WP**2 / (QM*N/L)
rho_back = -Q*N/L  # Background density?
dx = L / NG 

p = np.concatenate([np.arange(N), np.arange(N)])  # Particles? Note, the indices in MATLAB start from 1
Poisson = sparse.spdiags(([1, -2, 1] * np.ones((1, NG-1), dtype=int).T).T, [-1, 0, 1], NG-1, NG-1)
Poisson = Poisson.tocsc()

# electrons
xp = np.linspace(0, L-L/N, N).T   # initial positions of particles
pp = VT * (1 - VT**2)**(-0.5) * np.random.randn(N) # maxwellian   # Initial momentum
pm = np.arange(N)  # was transposed
pm = 1 - 2 * np.mod(pm+1, 2)
pp += pm * (V0 * (1 - V0**2)**(-0.5)) # Momentum + perturbation
gamma = (1. + pp**2)**0.5
vp = pp / gamma  # Velocity perturbation

# electron perturbation
xp = xp + XP1 * (L/N) * np.sin(2 * np.pi * xp / L * mode)
xp[np.where(xp < 0)] += L
xp[np.where(xp >= L)] -= L

histEnergy, histPotE, histKinE, histMomentum, t = [], [], [], [], []

if verbose:
    plt.figure(1)

for it in xrange(NT):   
    # Updating xp: momentum
    xp += vp * DT
    # Enforce periodicity
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
    #rho = Q / dx * mat.sum(axis=0) + rho_back 

    # Compute field
    #Phi = np.linalg.lstsq(Poisson.toarray(), -dx**2 * rho[:NG-1])[0]
    Phi = linalg.spsolve(Poisson, -dx**2 * rho[0:NG-1])
    Phi = np.concatenate((Phi,[0]))
    Eg = (np.roll(Phi, 1) - np.roll(Phi, -1)) / (2*dx)
  
    # Compute gamma'
    # Project q->p e aggiornamento velocita''
    pp = pp + mat*QM*Eg*DT
    gamma = (1. + pp**2)**0.5
    vp = pp / gamma
   
    Etot = Q/QM * (gamma-1).sum() + 0.5 * (Eg**2).sum() * dx
    histEnergy.append(Etot)
    histPotE.append(0.5 * (Eg**2).sum() * dx)
    histKinE.append(Q/QM * gamma.sum())
    histMomentum.append(Q/QM * pp.sum())
    t.append(it*DT)
   
    if (np.mod(it, TOut) == 0) and verbose:
        # Phase space        
        plt.clf()
        plt.subplot(2, 2, 1)
        plt.scatter(xp[0:-1:2], pp[0:-1:2], marker='.', color='blue')
        plt.scatter(xp[1:-1:2], pp[1:-1:2], marker='.', color='red')
        plt.xlabel('X')
        plt.ylabel('P')
        #plt.legend(loc=1)
        plt.title('$\omega_{pe}t=$' + str(DT*it))
	# Electric field
        plt.subplot(2, 2, 2)
        xg = np.linspace(0, L-dx, NG) + dx/2
        plt.plot(xg, Eg, label='E')
        plt.legend(loc=1)
	# Energies
        plt.subplot(2, 2, 3)
        plt.plot(t, histEnergy, label='Energy')
        plt.plot(t, histPotE, label='Potential')
        plt.plot(t, histKinE, label='Kinetic')
        plt.legend(loc=3)
        # Momentum
        plt.subplot(2, 2, 4)
        plt.plot(t, histMomentum, label='Momentum')
        plt.legend(loc=1)
        plt.pause(0.1)
        print it

print 'Time elapsed: ', time.clock() - start_time
#save duesmdt time histEnergy
if verbose:
    plt.figure(2)
    plt.plot(t, (histEnergy-histEnergy[0])/histEnergy[0])
    plt.xlabel('$\omega_{pe}t$')
    plt.ylabel('$(E(t)-E(0))/E(0)$')
    plt.pause(0.01)
    plt.show()
