# Rayleigh-Benard convection.

Parameters are set in routine rayleigh_const for convenience.

With this nondimensionalization, set rho==1 (parameter p1 in .rea
file) and visocity (p2) to be the desired Prandtl number.

Rayleigh number is set as Ra = Rc*(1+eps),  Rc=p76, eps=p75.

The buoyancy is ffy = Ra Pr T, where T is determined by 
boundary and initial conditions.

Critical Rayleigh number is around 1707.762 
(Somehow I was recalling 1734, but that appears to be for a 
particular geometric configuration considered by Laurette Tuckerman 
& Dwight Barkley)

GEOMETRY:

There are two primary cases, ray1.box and ray2.box.
The former specifies 10 elements in x, the latter only 9,
both for a 9x1 domain.

NOTES:

A time trace of (u,v)_max vs t is output to the logfile. See userchk.

Be careful about selecting an even number of elements in x
as it appears that the RB system likes to lock onto the grid spacing
and give a number of rolls that matches the number of elements, if the
elements have order-unity aspect ratio, as in the present case.
Thus, in the case, the 9 element mesh is likely to be more faithful
to the linear stability theory, at least for modest polynomial orders
of lx1=12.

It appears that one cannot realize Courant conditions of CFL ~ 0.5
with these cases because of the explicit convection treatment.
The given value dt=.02 is stable with lx1=12.

