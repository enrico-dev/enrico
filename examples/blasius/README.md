# Blasius boundary layer test case in 2D.

It tests the "ON " boundary condition (Outflow, Normal only), given by:

     u_parallel  - specified in userbc  (ux, uy, or uz)
     d (u_n)/ dn = 0

Note that the ON bc works only if the surface in question is normal
to either the X, Y, or Z axis.

In this test case the outflow condition is given by "O  ", which 
implies p + dU/dn = 0.

This configuration is set up for ymax=10, xmax ~= 480.   

At the end of the timestepping, most of the error is near the
outflow boundary, as would be expected.

The first field that is dumped is the exact solution, which is
also the initial condition. Thus, you can check the error by
plotting the difference between the first time slice and any 
later slice.