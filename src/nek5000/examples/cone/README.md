# Convecting cone problem

These cases provide examples of the convecting cone problem illustrated 
in Gottlieb & Orszg (1972) and Deville, Fischer, & Mund (2002).

Typically, to run them you would change SIZE such that the product
lx1*E1 = 32, where the total number of elements in the given case
is E=E1*E1.   For cone016 we have E=16 and E1=4, etc.  Thus, set lx1
as follows:

    cone256    lx1=2   
    cone064    lx1=4   
    cone016    lx1=8   

In this setup, each case (cone016, cone064, cone256) has its own
directory with its SIZE, .usr, and .box file.  To run any one of 
these examples, cd to the desired case and run genbox, to produce 
the .rea file, followed by genmap, to create the .map file). 

Note, for this problem, which is a pure convection problem, all of
the cost is in evaluation of the convection operator.  Thus, it 
does not make sense to use the characteristics scheme, which 
allows for a larger timestep at the expense of an increased 
number of evaluations of the convection operator.  (Generally,
the characteristics scheme is useful when the pressure solution
costs are high.)   Hence, in this example, we have IFCHAR = F.

