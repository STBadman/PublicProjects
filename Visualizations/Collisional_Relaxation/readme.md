max_relax.py

Relaxation to a Maxwellian Simulation

Samuel Badman - 12.03.17

--------------------------------------------------------------------------------

-------
Intro :
-------

This is a python script (tested in python 3) which simulates how 2 counter-
propagating beams of neutral particles which subsequently undergo elastic 
collisions relax to a maxwellian distribution of speeds. Users can modify 
the number of particles (default N=70) and the the radii of the particles
(default r = 3e-8 cm) by supplying these numbers as the first and second
command line arguments. The particle positions in the box are plotted in the
left hand panel, and the speed distribution is plotted in the right hand panel. 
Both an instantaneous (spiky) speed distribution, and a running average of the 
distributions are plotted.

--------------------------------------------------------------------------------

---------------------
Instructions to run :
---------------------

Step 1 - Download the python script max_relax.py to a local directory

Step 2 - Ensure you have python installed with numpy and matplotlib

Step 3 - Run the script from command line:


----------------------
Mac/Linux (Terminal) :
----------------------

(Default parameters)
$> cd /path/to/max_relax.py
$> python max_relax.py

(User supplied parameters) 
$> cd /path/to/script
$> python max_relax.py N r

N - number of particles in box : must be an integer, e.g. 100 
r - radius of particles        : must be a float expressed in cm, e.g. 2e-8
                                 N.B. the box dimensions are 1e-6cm so particle
                                 radius should be much smaller than this.
                                 
Both arguments should be supplied in the above order or defaults will be used.

--------------------------                                 
Windows (Command Prompt) :
--------------------------

C:\path\to\python.exe C:\path\to\max_relax.py                                

--------------------------------------------------------------------------------

---------                                 
Details :
---------

Box:
The simulation takes place in a 2D box of 10nm x 10nm dimensions with periodic
boundary conditions, i.e. a particle travelling off the right hand side of the 
box immediately appears on the left hand side with momentum conserved.

Collisions:
The particles are assumed to be hard spheres, and the collisions are modelled 
as causing the colliding particles momentum vectors to be reflected in the plane
tangent to the particles surfaces as viewed in the particles center of mass 
frame. 

Particle Properties:
All particles have the same mass - the mass of a hydrogen atom. 
The initial speeds of the particles are chosen by setting a  beam temperature 
of 300K then requiring v = sqrt(2kT/m). All particles have the same radius. All
particles are hard, neutral spheres which collide with each other elastically. 
50% of the particles are given initial speed +v and are initialized to lie in 
the left hand side of the box, while the other 50% are given initial speed -v
and initialized in the right hand side of the box.

Simulation:
The simulation is run with a timestep equal to the time it would take 1 
particle in the beam to cross 1/50th the width of the box. The simulation is 
run for 1000 timesteps.

Plotting:
At each timestep in the simulation, the results are plotted live on the screen,
producing an animation. There are two panels: (1) showing the positions of the 
the particle in the box. This allows the user to visually follow the behaviour
of the particles. The two beams are intially seen to move unimpeded but as 
collisions occur, the particle motion clearly isotropizes. (2) Showing the speed 
distribution which begins as a delta function at speed v, and relaxes over time
into the maxwellian shape which (in 2D) is proportional to |v|exp(-mv^2/2k_bT).
The instantaneous distribution at each timestep is plotted (this is spiky due
to random instantaneous fluctuations), and a running average which shows how the
initial delta function gradually smooths into the broad maxwellian distribution
we expect. The average is only started after 20 timesteps to stop the initially
very spiky distribution having too much weighting in the average at late times.
Finally, the analytic expectation of |v|exp(-mv^2/2k_bT) is plotted in red and
we can see how the averaged distribution converges to this curve. 
