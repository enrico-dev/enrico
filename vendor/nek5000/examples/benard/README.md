# Rayleigh-Benard convection.

Parameters are set in routine rayleigh_const for convenience.

With this nondimensionalization, set rho==1 (parameter p1 in .rea
file) and visocity (p2) to be the desired Prandtl number.

Rayleigh number is set as Ra = Rc*(1+eps),  Rc=p76, eps=p75.
Use p122 (steady state tolerance) sstol>0 to fit Ra_critical from
two additional Ra solutions prescribed by eps2 & eps3 (p121-122)

The buoyancy is ffy = Ra Pr T, where T is determined by 
boundary and initial conditions.

Critical Rayleigh number is around 1707.762 for the
wavenumber 3.117 [1].

NOTES:

A time trace of volume-avearged kinetic energy Ek vs t is output to
the logfile. See userchk and runtimeavg.

Be careful about selecting an even number of elements in x
as it appears that the RB system likes to lock onto the grid spacing
and give a number of rolls that matches the number of elements, if the
elements have order-unity aspect ratio, as in the present case.
Thus, in the case, the 9 element mesh is likely to be more faithful
to the linear stability theory, at least for modest polynomial orders
of lx1=8.

It appears that one cannot realize Courant conditions of CFL ~ 0.5
with these cases because of the explicit Boussinesq treatment.
The given value dt=.02 is stable with lx1=8.

Use 'grep Ra logfile' to get critical Ra or 'grep _ra logfile' to
get steady state amplitudes for each of 3 eps* cases in adition to
the last line w/ p75

[1] Chandrasekhar 1961 "Hydrodynamic and Hydromagnetic Stability", Table III on p. 43