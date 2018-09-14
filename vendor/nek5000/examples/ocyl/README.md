# Oscillating cylinder using ALE.


This case illustrates the Nek5000 ALE (arbitrary Lagrangian-Eulerian)
formulation for flow past an oscillating cylinder in 2D.

There are two cases:

* ocyl,  which uses Nek5000's linear elasticity solver to compute the mesh velocity at each step, and

* ocyl2, which uses a _user_supplied_ routine to generate an interpolating function for the ALE formulation.  

The key for the ALE is that the mesh velocity match the fluid velocity
at the domain boundaries (more precisely, that their normal component
match) and that the mesh velocity be as smooth as possible throughout
the domain.   Thus, determination of the mesh velocity simply requires
finding a smooth interpolant that will blend between any two parts of 
the domain boundary.   In ocyl2, this is done by solve a variable coefficient
Poisson problem that forces most of the blending to take place away from
the region of interest.  

In these examples, ocyl2 will be _much_ faster because there is no solve
for the mesh velocity on each step.

Frequency and amplitude of cylinder motion are set in the
.rea file:

  0.200000     p033  user: frequency of oscillation
  0.500000     p034  user: amplitude of oscillation

Note that the amplitude cannot be too large or the
mesh will become entangled and the simulation will stop.


The ALE formulation in this case solves an elasticity
problem at every timestep in order compute the mesh 
velocity.  While robust, this approach is slow and not
recommended for production runs.   A better approach is
to compute the mesh velocity for a single periodic cycle
and to then use a short Fourier exapansion based on a 
small set of time points over the full cycle to reconstruct
the desired mesh velocity at any future time point.
Such an approach would avoid repetitious computation of
the mesh velocity, which for this example is far more expensive
than computing the velocity/pressure for the fluid.  There
would be no loss in accuracy because the ALE formulation
allows for arbitrary mesh velocities, provided that they
meet the kinematic condition that the normal mesh velocity
component coincides with the normal fluid
velocity component on the boundary.


Prescribed boundary conditions are given in Cartesian coordinates.  
The prescribed surface motion is indicated by the "mv " bc on your 
moving surface (i.e., the cylinder surface). 

To use mesh motion, you must also have

T ifstrs
T ifmvbd 

in the .rea file and SIZE must have the following settings:

      parameter (lx2=lx1-2)
      parameter (ly2=ly1-2)
      parameter (lz2=lz1  )
      parameter (lx1m=lx1,ly1m=ly1,lz1m=lz1)
