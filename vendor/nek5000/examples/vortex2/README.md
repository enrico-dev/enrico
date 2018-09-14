# Vortex breakdown in a container with a rotating lid (axisymmetric case).


See vortex example for more detail about the corresponding 3D simulation 
configuration at Re=1854.

This is an axisymmetric example of the vortex breakdown problem
described in the Nek5000 Examples document.   

The domain size is 2x1.  The rotating lid is located at x=2.
For Re=1492, there is a single flow reversal (bubble) on the axis.

In this example, the azimuthal velocity is stored in the temperature
array t(...).

To have temperature serve the role of azimuthal velocity, you must
set the flag IFAZIV to be true either in the .usr file or in the .rea
file.

Note that the .rea and .box files specify "A  " on the y=0 boundary,
which is requisite for any axisymmetric flow problem that includes
the x-axis.

Also, any element touching the x-axis must be oriented with the "r"
direction pointing in the positive x direction.  (This condition
is automatically met for elements generated with genbox.)