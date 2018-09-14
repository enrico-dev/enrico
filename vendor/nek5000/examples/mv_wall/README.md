# 2D channel with a moving indentation.
  
This case demonstrates the Nek5000 ALE formulation for 
flow in a 2D channel with a moving indentation.  
 
This example has been the subject of several experimental 
and numerical studies [1,2,3].  
 
Illustration:

	                              Fixed Wall
	  ----------------------------------------------------------------
	                                     ^ 
	                                     |
	Inlet                                | b                            Outlet
	  -->                                |                                 
	               ____________          |                               
	              /            \         v
	  ___________/   moving     \_____________________________________
                   wall 
	 		^
	     y |
	 		|--->
	  		 x

A steady parabolic velocity profile is prescribed at the inlet. An 
indentation on the bottom wall moves in and out sinusoidally where
it's retracted position is flush with the wall.  
 
There are two main non-dimensional parameters given by: 
 
1) Reynolds Number, Re = U*b/nu = 507
2) Strouhal Number, St = b / (U*T) = 0.037

where,

	U := Mean Inlet Velocity
	b := The channel height
	T := The oscillation period

he motion of the bottom moving wall is given by the function

	y = F(x,t) = g(x)*h(t)    
 
where,

	           ----
	           | 
	           | 0 for t < 0
	  h(t)  =  |
	           | 0.5*[1 - cos(2*pi*t/T)]  for t < 0
	           |
	           ----


	g(x)  = 0.5*(h_max)*[1 - tanh (a |x| - x_2)]

With parameters from [1,2]

	h_max = 0.38*b
	a     = 4.14
	x_2   = 5.25*b

The entire length of the channel is given by, L = 40*b.  

This case is currently setup for Re = 507 and St = 0.037 in order to
compare with the experiments [1,2].
 
The ALE formulation for this problem proceeds by solving for
the mesh velocity via the Laplace equation, which
relies on the maximum principle to give a bounded interpolant. This ensures
that the boundary layer elements near the walls are preserved. The bulk of the
mesh deformation is pushed into the far field.   

This case also uses the turb_outflow boundary condition to ensure for 
characteristics to exit the domain at the outflow boundary.

[1] T.J. Pedley and K.D. Stephanoff, "Flow along a channel with 
a time-dependent indentation in one wall: the generation of
of vorticity waves," in Journal of Fluid Mechanics (1985), 
vol. 160, pp. 337-376.

[2] M.E. Ralph and T.J. Pedley, "Flow in a Channel with 
a moving indentation," in Journal of Fluid Mechanics (1988), 
vol. 190, pp. 87-112.
              _ ,            ,
[3] I. Demirdzic and M. Peric, "Finite volume method for prediction
of fluid flow in arbitrarily shaped domains with moving
boundaries," in International Journal for 
Numerical Methods in Fluids (1990), vol. 10, pp. 771-790.
