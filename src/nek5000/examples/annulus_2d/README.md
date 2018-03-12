# Natural convection between the two concentric cylinders

This case demonstrates the application of Nek5000 to simulation of natural
convection between the two concentric cylinders as considered by Grigull \& Hauf [1]
and presented by Van Dyke [2].

The inner cylinder (diameter _D/3_) is slightly heated with with 
respect to the outer one (diameter _D_). The Boussinesq approximation is used
to formulate the equations of motion, valid in situations where density
differences are small enough to be neglected everywhere except in the
gravitational forcing.  

Normalizing the Navier-Stokes and energy equations with 
* _D_ for the length scale 
* _D/U_ for the time scale (_U_ is the characteristic velocity in the given problem)
* _(T-T\_0)/(T_1-T\_0)_ nondimensional temperature, where _T\_0_ and _T\_1_ are the
respective temperatures of the outer and inner cylinders.
The main nondimensional parameters are:
* _Gr_ Grashof number
* _Pr_  Prandtl number

The simulations are performed for _Gr_ = 120000 and _Pr_=0.8 
until a steady-state solution is obtained, indicated by the change in velocity
magnitude between successive time steps becomes sufficiently small.

[1] Grigull \& Hauf, _Proc. of the 3rd Int. Heat Transfer Conf. 2_, p. 182--195 (1966)

[2] M. Van Dyke _An Album of Fluid Motion_, Parabolic Press, Stanford, CA, (1982)