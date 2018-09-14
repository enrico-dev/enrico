# Double shear-layer rol-up

Thin and thick shear layer cases are considered.

Also, these cases illustrate transport of multiple passive scalars.

This is the .usr file for the double shear-layer roll-up problem
that has been studied by many authors, including

    Bell, Collela, & Glaz, JCP [85],  p257 (1989)
    Brown and Minion,      JCP [122], p165 (1995)
    Tadmor                 1997, Italy
    C.W. Shu               p.65, 1997 ICASE Rep.
    Fischer & Mullen       CRAS (2001)

The challenge presented by this problem is that the shear
layers thin out with time, eventually evolving to a state
where they are thinner than the grid resolution.  Ideally,
a code should be able to handle this situation with "graceful
degradation" in which the macroscopic properties continue to
evolve without blow-up of the solution.

Fischer & Mullen (CRAS 2001) illustrate that the filter-based
stabilization is an effective approach to solving this problem.
Further investigation pointed to aliasing of the nonlinear
quadrature as a major stability issue in the original SEM
formulations.   By over-integrating the convective transport
term one can recover stability for this problem and run without
filtering and even with arbitrarily small viscosity nu > 0.
Thus, in the .rea file, p99=3 (dealiasing on) and p103=0
(filtering off).

There are two different initial conditions - the "thin" shear
layer (scaling parameter rho=100) and the thick one (rho=30),
with the thin case considerably harder than the thick.  
Reasonable results for the thin case can be realized 
with 256^2 resolution (lx1=17 with the present 16x16 element mesh).
For results with the thin-layer case, see F&M (2001).
Reasonable results for the thick case can be realized 
with 128^2 resolution (lx1=9 with the present 16x16 element mesh).

Reynolds numbers of 10^4 and 10^5 have been reported in the
literature.  In the present case, the .rea file has param 2
set to -100000, corresponding to nu=.00001 and Re=10^5.

# PASSIVE SCALAR TRANSPORT

The .rea and SIZEu files are also configured for the transport
of 4 passive scalars, just to illustrate how this is accomplished.
The scalar fields are initialized with some arbitrary values in
useric.  The conductivities are set to 1, .1, .01, and .001,
with the first one set as the thermal conductivity (param 8)
and the others set on line 124 in the .rea file:

	1.e-1         1.e-2         1.e-3         1.00000       1.00000    

Note that Nek5000 supports, in principle, an arbitrary number
of passive scalars (up to 20 have been used to date).   Use 
the .usr file to initialize the properties for a very high 
number of passive scalars.  In addition, one must increase 
the "ldimt" parameter in the SIZE file and recompile the 
entire code.
