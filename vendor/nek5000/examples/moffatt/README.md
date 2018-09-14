# Moffatt eddies - Viscous Flow in a Wedge.

Important files here:

`SIZE`

`moffatt`   case files

`moff_circ` case files

`moff2.m`

`run_all`   -- a short script to build and run the tests.

## moffatt

In [1], Moffatt gives the asymptotic size and strength of vortices
("Moffatt eddies") for Stokes flow in wedge shaped domains that were 
later visualized experimentally by Taneda [2].

Equations (3.11a) and (3.12) in [1] give the respective ratio
of successive eddy sizes and strengths, as a function of parameters
p1 and q1, satisfying

	p1 := xi  / (2 alpha)
	q1 := eta / (2 alpha)

	sin xi  cosh eta  + k xi    = 0           (3.6)
	cos xi  sinh eta  + k eta   = 0

	k := sin (2 alpha) / (2 alpha),

where 2 alpha is the included interior angle of the wedge domain.

Aside from xi=eta=0, (3.6) has solutions near (xi,eta)=(4.21,2.26),
as indicated in Table 1 of [1].

Using these values as starting points, we find that the pair

	(xi,eta)=(4.229163797472673 ,  2.160973066794531 )

is a valid solution for 2*alpha=28.5 degrees, as generated in the
matlab script, "moff2.m"

The corresponding asymptotic size and strength ratios are:

	rn/r_n+1 = exp(pi/q1):      size_ratio = 2.030997533604236e+00 .
	vn/v_n+1 = exp(pi*p1/q1):   strg_ratio = 4.083778044216303e+02 .

Using Nek5000's steady Stokes solver (based on the Uzawa algorithm
described in [3]) this "moff" case generates several vortices whose
size and strength is in excellent agreemeent with the results of [1].

The vortices are quite weak near the wedge vertex, but the strongest
ones near the top, where there is sufficient signal to noise, agree
well with the respective strength and size ratios, 408.378 and 2.031.

The results in the table below are produced by this usr file (followed
by "grep eddy logfile").

Eddy  Crossing point  peak velocity    size ratio     strength ratio

	14   2.7600000E-02   8.2902315E-12   3.6315789E+00   3.3004314E+02 eddy
	15   5.5200000E-02   3.3646379E-09   2.0000000E+00   4.0585572E+02 eddy
	16   1.1160000E-01   1.3740772E-06   2.0217391E+00   4.0838783E+02 eddy
	17   2.2720000E-01   5.6113675E-04   2.0358423E+00   4.0837354E+02 eddy
	18   4.5240000E-01   2.1748396E-01   1.9911972E+00   3.8757747E+02 eddy

Because the strenght of the eddies falls off so quickly, it's imperative
to take steps to minimize noise when generating and post-processing the solution.
In this case, relatively few elements were used, but the order is N=19. The
Uzawa scheme is based on the Pn-Pn-2 formulation, so we have N=17 for pressure,
meaning that (lx1=20, lx2=lx1-2).

The iteration tolerances are very tight: p20=1e-20, p21=1e-21 for pressure
and velocity.


## moff_circ

This is a Moffatt-eddy case [1] with a circular lid, similar to [2].


## Nek5000 test mode:

If this case fails (i.e., does not yield a sufficiently accurate answer), it will print:

      ERROR TOO LARGE FOR MOFFATT CASE   err, tol

otherwise it will print:

       Success: relative strength/size errors:   8.6596E-06  2.3842E-03

[1]  "Viscous and resistive eddies near a sharp corner" H.K. Moffatt, J. Fluid Mech., 
(18) 1, January, 1964, pp. 1--18.

[2]  "Visualization of Separating Stokes Flows" Sadatoshi Taneda, J. Phys. Soc. Jpn., 
(46) 6, June, 1979 pp. 1935--1942.

[3]   R{\o}nquist, E., "Optimal Spectral Element Methods for the Unsteady
Three-Dimensional Incompressible {Navier-Stokes} Equations", 
PHD thesis, Massachusetts Institute of Technology, 1988

