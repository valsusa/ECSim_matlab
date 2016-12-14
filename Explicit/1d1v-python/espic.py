# -*- coding: utf-8 -*-

""" ESPIC: a simple 1D1V electrostatic PIC code 
    Based on 'espic.c' developed for 2000 Heraeus School, Jena
  
  (c) P. Gibbon, KU-Leuven & FZ-Juelich, November 2013

"""
import numpy as np  ## numeric routines; arrays
import pylab as plt  ## plotting

#include "es.h"			/* common variables 

""" core routines
"""

#==============================================================================
#  Particle loading - positions 
#==============================================================================
    
def loadx(bc_particle):            
    global dx, grid_length, rho0, npart, q_over_me, a0
    global charge,mass,wall_left,wall_right
    print "Load particles"

# set up particle limits
    if (bc_particle >= 2):  
        # reflective boundaries
        # place particle walls half a mesh spacing inside field boundaries 
        wall_left = dx/2.
        wall_right = grid_length-3*dx/2.
        plasma_start = wall_left
        plasma_end = wall_right   # actually want min(end,wr) */
     
    else:
        # periodic boundaries
        plasma_start = 0.
        plasma_end = grid_length
        wall_left = 0.
        wall_right = grid_length
      
    xload = plasma_end - plasma_start  # length for particle loading */
    dpx = xload/npart                    #  particle spacing */
    charge = -rho0*dpx         #  pseudo-particle charge normalised to give ncrit=1 (rhoc=-1) 
    mass = charge/q_over_me  #  pseudo-particle mass (need for kinetic energy diagnostic) 

    for i in range(npart):
        x[i] = plasma_start + dpx*(i+0.5)  # Python ndarrays start at index 0
 #       x[i] += a0*np.cos(x[i]) # Include small perturbation
        
    return True

    
    
#==============================================================================
#  Particle loading - velocities
#==============================================================================
    
def loadv(idist,vte):       
    global npart,v
    print "Set up velocity distribution"
    iseed = 76523  # random number seeds 
    idum1 = 137 

    if ( idist == 1 ):
       # >10 = set up thermal distribution   
        for i in range(npart):
            vm = vte*np.sqrt( (-2.*np.log((i+0.5)/npart)) ) # inverted 2v-distribution - amplitude */
            rs = np.random.random_sample() # random angle */
            theta = 2*np.pi*rs
            v[i] = vm*np.sin(theta) # x-component of v 
   
      # scramble particle indicies to remove correlations between x and v
        np.random.shuffle(v)

    else:  
        # Default is cold plasma */
         v[1:npart] = 0.

    return True
    
#==============================================================================
#   Compute densities
#==============================================================================

def density(bc_field,qe):      
   global x,rhoe,rhoi,dx,npart,ngrid,wall_left,wall_right
   j1=np.dtype(np.int32) 
   j2=np.dtype(np.int32) 

   re = qe/dx #  charge weighting factor 
   rhoe=np.zeros(ngrid+1)        # electron density 
    # map charges onto grid
   for i in range(npart):
      xa = x[i]/dx
      j1 = int(xa)
      j2 = j1 + 1
      f2 = xa - j1
      f1 = 1.0 - f2
      rhoe[j1] = rhoe[j1] + re*f1
      rhoe[j2] = rhoe[j2] + re*f2

   if (bc_field == 1):
      #  periodic boundaries */ 
      rhoe[0] += rhoe[ngrid]
      rhoe[ngrid] = rhoe[0]

   elif (bc_field == 2):
     #  reflective - 1st and last (ghost) cells folded back onto physical grid */
      iwl = wall_left/dx
      rhoe[iwl+1] += rhoe[iwl]
      rhoe[iwl] = 0.0
      iwr = wall_right/dx
      rhoe[iwr] += rhoe[iwr+1]
      rhoe[iwr+1] = rhoe[iwr]
   else:
      print "Invalid value for bc_field:", bc_field
      
#  Add neutral ion density      
   rhoi = rho0
#   print rhoe[0:ngrid+1]
   return True

#==============================================================================
#  Compute electrostatic field
#==============================================================================

def field():          
  global rhoe,rhoi,ex,dx,ngrid
  
  rhot=rhoe+rhoi  # net charge density on grid
  
  # integrate div.E=rho directly (trapezium approx) 
  # end point - ex=0 mirror at right wall 

  Ex[ngrid]=0.          # zero electric field  
  edc = 0.0

  for j in range(ngrid-1,-1,-1):
      Ex[j] = Ex[j+1] - 0.5*( rhot[j] + rhot[j+1] )*dx
      edc = edc + Ex[j]

  if (bc_field == 1):
	# periodic fields:  subtract off DC component */ 
	#  -- need this for consistency with charge conservation */
      Ex[0:ngrid] -= edc/ngrid  
      Ex[ngrid] = Ex[0]
      
  return True
    
    
#==============================================================================
#  Particle pusher
#==============================================================================

def push():  
    global x,v,Ex,dt,dx,npart,q_over_me

    for i in range(npart):

      #  interpolate field Ex from grid to particle */
      xa = x[i]/dx
      j1 = int(xa)
      j2 = j1 + 1
      b2 = xa - j1
      b1 = 1.0 - b2
      exi = b1*Ex[j1] + b2*Ex[j2]
      v[i] = v[i] + q_over_me*dt*exi       #  update velocities */

    
    x += dt*v    #  update positions (2nd half of leap-frog) 
    
    return True
    
    
#==============================================================================
#  check particle boundary conditions
#==============================================================================
    
def particle_bc(bc_particle,xl):      
   global x
#  int iseed1 = 28631;        /* random number seed */
#  int iseed2 = 1631;         /* random number seed */
 

  #  loop over all particles to see if any have
  #     left simulation region: if so, we put them back again
  #     according to the switch 'bc_particle' **/

   if (bc_particle == 1):
     #  periodic
     for i in range(npart):
       if ( x[i] < 0.0 ):
         x[i] += xl
       elif ( x[i] >= xl ):
         x[i] -= xl

   elif (bc_particle == 2):
     # reflective
     print "Reflective boundaries not implemented yet."

   elif (bc_particle == 3):
     # reflective
     print "Thermal boundaries not implemented yet."

   else:
     print "Invalid value for bc_particle:", bc_particle

   return True
    
    
#==============================================================================
#  Diagnostic outputs for fields and particles
#==============================================================================
    
def diagnostics():      
    global rhoe,Ex,ngrid,itime,grid_length,rho0,a0
    global ukin, upot, utot, udrift, utherm, emax,fv,fm
    global iout,igraph,iphase,ivdist
    xgrid=dx*np.arange(ngrid+1)
    if (itime==0): 
        plt.figure('fields')
        plt.clf()
    if (igraph > 0):
      if (np.fmod(itime,igraph)==0): # plots every igraph steps
    # Net density
        plt.subplot(2, 2, 1)
        if (itime >0 ): plt.cla()
        plt.plot(xgrid, rhoe+rho0, 'r', label='density')
        plt.xlabel('x')
        plt.xlim(0,grid_length)
 #       plt.ylim(-2*a0,2*a0)
        plt.legend(loc=1)
    # Electric field
        plt.subplot(2, 2, 2)
        if (itime >0 ): plt.cla()
        plt.plot(xgrid, Ex, 'b', label='Ex')
        plt.xlabel('x')
#        plt.ylim(-2*a0,2*a0)
        plt.xlim(0,grid_length)
        plt.legend(loc=1)
        
        if (iphase > 0):
          if (np.fmod(itime,iphase)==0):  
    # Phase space plots every iphase steps
            axScatter = plt.subplot(2, 2, 3)
            if (itime >0 ): plt.cla()
            axScatter.scatter(x,v,marker='.',s=1)
    #    axScatter.set_aspect(1.)
            axScatter.set_xlim(0,grid_length)
            axScatter.set_ylim(-vmax,vmax)
            axScatter.set_xlabel('x')
            axScatter.set_ylabel('v')

        if (ivdist > 0):
          if (np.fmod(itime,ivdist)==0):
    # Distribution function plots every ivdist steps
            fv=np.zeros(nvbin+1)        # zero distn fn 
            dv = 2*vmax/nvbin   # bin separation */
            for i in range(npart):
                vax= ( v[i] + vmax )/dv   # norm. velocity */
                iv = int(vax)+1 # bin index */	 
                if (iv <= nvbin and iv > 0): fv[iv] +=1 #  /* increment dist. fn if within range
        
            plt.subplot(2, 2, 4)
            if (itime >0 ): plt.cla()
            vgrid=dv*np.arange(nvbin+1)-vmax
            plt.plot(vgrid, fv, 'g', label='f(v)')
            plt.xlabel('v')
            plt.xlim(-vmax,vmax)
#        plt.ylim(-2*a0,2*a0)
            plt.legend(loc=1)

    
        plt.pause(0.0001)
        plt.draw()
        filename = 'fields_%0*d'%(5, itime)
        if (iout > 0):
          if (np.fmod(itime,iout)==0):  # printed plots every iout steps
            plt.savefig(filename+'.png')

#   total kinetic energy
    v2=v**2
    vdrift=sum(v)/npart
    ukin[itime] = 0.5*mass*sum(v2)
    udrift[itime] = 0.5*mass*vdrift*vdrift*npart  
    utherm[itime] = ukin[itime] - udrift[itime]
  
#   potential energy 
    e2=Ex**2
    upot[itime] = 0.5*dx*sum(e2)
    emax = max(Ex) # max field for instability analysis */
 
#  total energy 
    utot[itime] = upot[itime] + ukin[itime]
    
    return True
  
#==============================================================================
#    Plot time-histories
#==============================================================================
  
def histories():

#FILE *history_file;     /* file for writing out time histories */            

  global ukin, upot, utot, udrift, utherm
  xgrid=dt*np.arange(nsteps+1)
#  plt.clf()
  plt.figure('Energies')
#  plt.subplot(2, 2, 1)
  plt.plot(xgrid, upot, 'b', label='Upot')
  plt.plot(xgrid, ukin, 'r', label='Ukin')
  plt.plot(xgrid, utot, 'black', label='Utot')
#  plt.plot(xgrid, udrift, 'g', label='Udrift')
  plt.xlabel('t')
  plt.ylabel('Energy')

#  plt.xlim(0,grid_length)
#  plt.ylim(-2*a0,2*a0)
  plt.legend(loc=1)
  plt.savefig('energies.png')

#  /*  write energies out to file */
#   if (itime==1) {fprintf(history_file, " t, U_drift, U_therm, U_field, U_total, Emax\n");}
#   fprintf( history_file, "%f  %e %e  %e  %e %e\n", itime*dt, udrift, utherm, upot, utot, emax );



#==============================================================================
#  Main program
#==============================================================================
  

npart=10000           # particles
ngrid=100             # grid points 
nsteps=50             # timesteps

# particle arrays
x = np.zeros(npart)  # positions
v = np.zeros(npart) # velocities#    particle_bc() 		 

# grid arrays
rhoe=np.zeros(ngrid+1)        # electron density 
rhoi=np.zeros(ngrid+1)        # ion density  
Ex=np.zeros(ngrid+1)          # electric field  
phi=np.zeros(ngrid+1)         # potential 
# time histories
ukin=np.zeros(nsteps+1)
upot=np.zeros(nsteps+1)
utherm=np.zeros(nsteps+1)
udrift=np.zeros(nsteps+1)
utot=np.zeros(nsteps+1)
 
# Define main variables and defaults
# ----------------------------------
 
#ni = 0;          /* # ions (fixed) */
#grid_length = 4.0*np.pi  # size of spatial grid
grid_length = 10.0  # size of spatial grid
plasma_start = 0.           # LH plasma edge 
plasma_end = grid_length    # RH plasma edge
dx = grid_length/ngrid
dt = 0.1             # normalised timestep
q_over_me=-1.0       # electron charge:mass ratio
rho0 = 1.0           # background ion density 
vte = 0.001           # thermal velocity 
vmax = 5*vte       # max velocity for f(v) plot
nvbin=50            # bins for f(v) plot
a0 = 0.1             # perturbation amplitude
wall_left=0.
wall_right=1.
bc_field = 1       #  field boundary conditions:  1 = periodic
                    #                              2 = reflective 
bc_particle = 1     # particle BCs:  1 = periodic
                    #                2 = reflective
                    #                3 = thermal
profile = 1   	 # density profile switch 
distribution = 1  # velocity distribution switch: 0 = cold, 1=thermal
			 # 				          2 = 2-stream 
ihist = 5         # frequency of time-history output 
igraph = 1       # freq. of graphical snapshots 
iphase = 0
ivdist = 1
iout = 0        # freq. of saved graphics files 
itime = 0         # initialise time counter 

   
   
#  Setup initial particle distribution and fields
#  ----------------------------------------------
   
loadx(bc_particle)              # load particles onto grid 
loadv(distribution,vte)			# define velocity distribution 
x += 0.5*dt*v                   #  centre positions for 1st leap-frog step 
particle_bc(bc_particle,grid_length)
density(bc_field,charge)	  	 # compute initial density from particles 
field()			            # compute initial electric field 
diagnostics()		            # output initial conditions 
print 'resolution dx/\lambda_D=',dx/vte

#  Main iteration loop
#  -------------------

for itime in range(1,nsteps+1):
    print 'timestep ',itime
    push()			       # Push particles 
    particle_bc(bc_particle,grid_length) # enforce particle boundary conditions
    density(bc_field,charge)	    # compute density 
    field()		               # compute electric field (Poisson) 
    diagnostics()	# output snapshots and time-histories

histories()   # Produce time-history plots
      
#raw_input("Press key to finish")
#plt.close()
print 'done'
