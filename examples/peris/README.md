# Peristaltic flow.

This .usr/.rea/SIZE trio provides a rudimentary example of peristaltic
flow in a tube with inflow/outflow boundary conditions.  

Technically, this example is flawed because the amplitude and frequency
of the peristaltic pumping has no impact on the overall flow rate, due
to the Dirichlet inflow condition specified for velocity.  It thus
should be viewed only as an illustration of one approach to implementing
the boundary conditions and the mesh velocity for the ALE formulation.

This case also illustrates the use of the special outflow treatment,
where div U > 0 is imposed in the last layer of elements in order
to ensure that all characteristics are outgoing where outflow 
conditions have been specified.

The statement ifusermv=.true. in userdat2() indicates that the user
will provide the mesh velocity, wx,wy, and wz.   Note that in SIZE,
it is necessary to set lx1m=lx1, ly1m=ly1, and lz1m=lz1 to allocate
the moving mesh array space.