# Eddy solutions with moving boundary conditions.

This case monitors the error for an exact 2D solution
to the Navier-Stokes equations based on the paper of Walsh [1],
with an additional translational velocity (u0,v0) and a moving
domain.

The initial computational domain is [0,2pi]^2 with doubly-periodic
boundary conditions.

Walsh's solution consists of an array of vortices determined 
as a linear combinations of eigenfunctions of having form:

    cos(pi m x)cos(pi n y), cos(pi m x)sin(pi n y)
    sin(pi m x)cos(pi n y), sin(pi m x)sin(pi n y)

and

    cos(pi k x)cos(pi l y), cos(pi k x)sin(pi l y)
    sin(pi k x)cos(pi l y), sin(pi k x)sin(pi l y)

While there are constraints on admissible (m,n),(k,l)
pairings, Walsh shows that there is a large class of
possible pairings that give rise to very complex vortex
patterns.

Walsh's solution applies either to unsteady Stokes or 
unsteady Navier-Stokes.  The solution is a non-translating
decaying array of vortices that decays at the rate 

     exp ( -4 pi^2 (m^2+n^2) visc time ),

with (m^2+n^2) = (k^2+l^2). A nearly stationary state may
be obtained by taking the viscosity to be extremely small,
so the effective decay is negligible.   This limit, however,
leads to an unstable state, thus diminsishing the value of 
Walsh's solution as a high-Reynolds number test case.

It is possible to extend Walsh's solution to a stable convectively-
dominated case by simulating an array of vortices that translate
at arbitrary speed by adding a constant to the initial velocity field.  
This approach provides a good test for convection-diffusion dynamics.
In the current file set the translational velocity is specified as:

    U0 =[u0,v0] := [param(96),param(97)]    ( in the .rea file )

The approach can also be extended to incompressible MHD with non-unit
magnetic Prandtl number Pm.

[1] Owen Walsh, "Eddy Solutions of the Navier-Stokes Equations,"
in The Navier-Stokes Equations II - Theory and Numerical Methods,
Proceedings, Oberwolfach 1991, J.G. Heywood, K. Masuda,
R. Rautmann,  S.A. Solonnikov, Eds., Springer-Verlag, pp. 306--309 (1992).
