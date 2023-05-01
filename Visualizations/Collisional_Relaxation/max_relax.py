'''
max_relax.py
Python script which runs an onscreen simulation of two counter propagatsing 
beams of neutral particles in 2D which collide and via many elastic collisions,
relax their speed distribution to a Maxwellian. The simulation takes place in
a box of dimensions 10 nm x 10 nm, with periodic boundary conditions. User can 
choose to specify the number of particles and their radius in cm as the 1st and
2nd arguments of the script.

By Samuel Badman for Physics/Astronomy C207 Class Project, UC Berkeley, 12.03.17
'''


#Import required modules. Need numpy and matplotlib
from math import pi
import numpy as np
import numpy.linalg as npl
import matplotlib.pylab as plt
import sys, warnings
#Clean output for graphics
warnings.filterwarnings("ignore",".*GUI is implemented.*")

## User can provide number and radius of particles as command line arguments
args = sys.argv
# Default values for case where no command line arguments provided
if len(args) != 3 :
    # Particle Size
    r_0 = 3.e-8 #cm
    # Number of Particles
    N = 70
    print('Using default arguments: ')
    print('To supply user input, use: python max_relax.py <<num of particles'
          +'>> <<particle radius (in cm)>>')
# Otherwise, take user provided command line arguments    
else : 
    N,r_0 = int(args[1]),float(args[2])
    print('Using user input arguments: ')
    
print('Arg 1 : Number of Particles : '+str(N))
print('Arg 2 : Radius of particles : '+str(r_0)+' cm')    
print('Use CTRL-C to quit simulation')

#Fundamental Constants in CGS
e   = 4.80320427e-10# esu electron charge
m_e = 9.10938215e-28# g electron mass
m_H = 1.6737236e-24 # g hydrogen mass
k_B = 1.38064852e-16 # erg/K boltzmann constant
c   = 2.99792458e+10# cm/s speed of light

#Particle Mass - choose hydrogen mass
m = m_H

#Initial velocity distribution is chosen to be two counterpropagating beams of
#hydrogen atoms, each with a beam temperature of 300K
#Initial Beam Temperature
T_i = 300. #K
# Beam particle velocities (Factor of 2 from equipartition theorem)
v_i = (2.*k_B*T_i/m)**0.5

#Box Dimensions - simulation occurs in 2D in a box of dimensions given below. 
# The box has periodic boundary conditions so particles leaving one side enter
# on the other. This prevents loss of particles.
x_box = 1.e-6 #cm 
y_box = 1.e-6 #cm

#Position and Velocity Vectors for Each Particle
pos = np.zeros([2,N])
vel = np.zeros([2,N])

#Initialize 
#Randomize Initial Positions within box
pos[0,:int(N/2)] = np.random.uniform(0.,x_box/2.,int(N/2))
pos[0,int(N/2):] = np.random.uniform(x_box/2.,x_box,int(N/2))
pos[1,:] = np.random.uniform(0.,y_box,N)
#Start every particle in 2 counter propagating beams in the x direction.
vel[0,:int(N/2)] = v_i
vel[0,int(N/2):N] = -v_i
vel[1,:] = 0.
    
# Choose time resolution such that the initial beam would take 50 time steps to
# cross the x dimension of the box.  
dt = x_box/v_i/50. #s
# Simulate for 1000 timesteps
t_stop = 1000.*dt #s

#Create time array to iterate over
time_array = range(int(t_stop/dt))

#Initialize figure
plt.figure(figsize = (16,8))
plt.subplot(121)
plt.xlim([0.,x_box/1.e-7])
plt.ylim([0.,y_box/1.e-7]) 

hist_av_x = 0.
hist_av_y = 0.
hist_av_sp = 0.

for i in time_array :
    
    #update time
    t = float(i)*dt
    
    #update velocities (collision operator)
    #iterate through particle relative distances, if 2 particles intersect, 
    #apply impulse such that the particles momenta are reflected at the boundary
    #of intersection (elastic collision)
    for particle1 in range(N) :
        for particle2 in range(particle1+1,N) :
            #Get center-center separation vector, magnitude and direction
            sep_vec = pos[:,particle2]-pos[:,particle1]
            mag   = npl.norm(sep_vec)
            direc = sep_vec/mag
            #Get relative velocity of two particles
            vel_diff = vel[:,particle2]-vel[:,particle1]
            
            #Apply impulse if particles intersect, and their relative velocity
            #is bring the centers closer together. (This last condition ensures
            #particles escape from each other rather than colliding multiple 
            #times if they don't escape intesection in 1 timestep).
            if mag <= 2.*r_0 and np.dot(direc, vel_diff) <= 0. : 
                impulse = np.dot(vel[:,particle2]-vel[:,particle1],direc)*direc
                vel[:,particle1] = vel[:,particle1] + impulse  
                vel[:,particle2] = vel[:,particle2] - impulse
                
    # Update positions according to new velocities. 
    pos[0,:] = (pos[0,:] + vel[0,:]*dt) % x_box
    pos[1,:] = (pos[1,:] + vel[1,:]*dt) % y_box
    
    
    #Create speed histogram
    bins_sp = np.array(list(range(0,60)))/16.*v_i 
    hist_sp,bins_sp = np.histogram((vel[0,:]**2 + vel[1,:]**2)**0.5,bins=
                                    bins_sp)
    bin_sp_centers = 0.5*(bins_sp[1:]+bins_sp[:-1])
    
    #Update running average of histograms 
    if i >= 20.: hist_av_sp = (hist_av_sp*float(i) + hist_sp)/float(i+1)
    else : hist_av_sp = hist_sp
    
    #Get analytic maxwellian speed distribution and normalize:
    maxwellian = bin_sp_centers*np.exp(-0.5*m*bin_sp_centers**2/(k_B*T_i))
    maxwellian = N*maxwellian/np.sum(maxwellian)

    # Update Plot
    # Position Plot
    plt.title('Particle Positions, Time = '+str(int(t*1e12))+' ps')
    plt.xlabel('x / nm')
    plt.ylabel('y / nm')
    plt.xlim([0.,x_box/1.e-7])
    plt.ylim([0.,y_box/1.e-7])
	
    for particle in range(N) :
        plt.gca().add_patch(plt.Circle((pos[0,particle]/1.e-7,
                                        pos[1,particle]/1.e-7),
                                        radius=r_0/1.e-7))
                                            
    #Plot Velocity histograms,the running average and the analytic 2D maxwellian
    # for T = T_i
    plt.subplot(122)
    plt.title('Speed distribution')
    plt.xlabel('Speed / units of initial beam speed')
    plt.ylabel('Number of particles per unit speed')
    plt.ylim([0.,max(N/8.,float(max(hist_sp)))])
    plt.plot(bin_sp_centers/v_i,hist_av_sp,label = 'Averaged')
    plt.plot(bin_sp_centers/v_i,hist_sp,linestyle=':',color='green', 
    label = 'Instantaneous')
    plt.plot(bin_sp_centers/v_i,maxwellian, color='red',label='Analytic')
    plt.legend(loc='upper right')
    
    #Update plot then clear for next timestep.        
    plt.pause(0.000001)
    plt.cla()
    plt.subplot(121)
    plt.cla()
