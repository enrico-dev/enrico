# Stratified 2D flow.

These are Boussinesq-based stratified-flow examples that illustrate 
the blocking phenomena in which fluid particles resist change in 
elevation because of the potential energy change required to overcome 
strong stratification.  The examples are based on a discussion found 
in Tritton, "Physical Fluid Dynamics"  (Oxford University Press, 1988).
Results for these cases can be found in the Nek5000 "users.pdf" examples
file.

Two cases are considered iwth Prandtl number equal 1 and 1000.

In the simulations here, the hydrostatic pressure distribution is altered 
so that we solve for a field that is uniform (i.e., not varying in the 
y direction) at outflow, as this is the standard outflow condition in 
Nek5000.  To illustrate the point, consider that the Boussinesq approximation 
with a constant temperature field would give rise to a uniform body force in 
the negative y-direction, and that the resulting hydrostatic pressure field 
would vary linearly in y.   In the (linear) stratified flow case the forcing 
is of the form
 
       ffy = temp / Fr2
 
where Fr2 is the Froude number squared and the temperature field has for 
initial- and boundary-conditions a linear variation in y.   As it stands, 
such a forcing would lead to a hydrostatice pressure field that varies 
quadratically in y, and it is clear that p=0 at outflow would be incorrect.   
To account for the vertical loading, we modify the forcing term to be
 
       ffy = (temp-y) / Fr2
 
which implies that the computed pressure differs from the actual pressure 
by -.5*y*y/Fr2.  For the incompressible Navier-Stokes equations and the 
given boundary conditions this modification has no impact on the dynamics 
of the flow.
