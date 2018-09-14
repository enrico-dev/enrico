# Thick-wall hollow sphere under an internal pressure.

This is an example of using the linear elasticity steady 3D solid solver 
in Nek5000
 
The test case corresponds to a thick-wall hollow sphere under an internal 
pressure P (and zero external pressure) which has the analytical solution 
 
	u_r=P*r1^3*[(1-2nu)r+(1+nu)r2^3/(2r^2)]/E/(r2^3-r1^3)
	u_theta=u_phi=0
(see, for example, Fung "Foundations of Solid Mechanics")    

E is Young's modulus, nu is Poisson's ratio

	r1=0.5-inner radius, r2=1-outer radius 

Traction boundary conditions are used in the current example: 
normal stress = -P at the inner sphere (b.c. option 'S'), 
stress-free at the outer sphere (b.c. option 'O')
IFSTRS=T is necessary to use traction boundary conditions 

Displacement (Dirichlet) boundary conditions can be used as well
(corresponding boundary values for displacements are also in userbc)
To use - change 'S' and 'O' to 'v' in .rea file 

Standard NEK arrays for velocities (vx,vy,vz) in fluid formulation
contain displacements in solid formulation

Description of the SEM method for linear elasticity and convergence 
results for the current test case, as well as for the other cases, 
can be found at http://www.mcs.anl.gov/~peet/solid.pdf

All additional subroutines necessary for using linear elasticity 
solver are currently in .usr file and include
- steady_elast_solve, elast, elast2d_e, elast3d_e (elasticity routines)
- cg, solveM, axcg (conjugate gradient routines)
- getstress (traction boundary conditions routine)
- energy_norm (routine for calculating error in the energy norm)
