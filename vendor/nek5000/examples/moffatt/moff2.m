%
%    matlab script to solve for quantities related to Moffatt eddies,
%    described in 
%
%      "Viscous and resistive eddies near a sharp corner"
%       H.K. Moffatt, J. Fluid Mech., (18) 1, January, 1964, 
%       pp. 1--18.
%
%    See the Nek5000 user file, moff.usr, for more detail.
%
%
%    sin a2p cosh a2q + p sin a2 = 0   a2=28.5 deg 
%    cos a2p sinh a2q + q sin a2 = 0
%
%    k = sin a2/a2
%
%    sin x cosh y + kx = 0
%    cos x sinh y + ky = 0
%
%    (2n-1) pi < x_n < (2n-1/2) pi
%
format compact; format longe; hold off

x=2*ones(2,1); f=0*x;

J=eye(2);

a2 = pi*(28.5)/180;      %% 28.5 degree interior angle
k  = sin(a2)/a2;

xi = 4.2; eta = 2.2;  %% Initial guess from Moffatt 64, Table 1.

f(1) = sin(xi)*cosh(eta) + xi*k;
f(2) = cos(xi)*sinh(eta) + eta*k;

x=[xi;eta];

for iter=1:20;  % Solve nonlinear system with Broyden's method

   s = -J\f; 
   x = x+s;

   fo = f; xi=x(1);eta=x(2); 

   f(1) = sin(xi)*cosh(eta) + xi*k;
   f(2) = cos(xi)*sinh(eta) + eta*k;
%  [ x(1) x(2) f(1) f(2) ]

   y = f-fo;
   J = J - ((J*s-y)*s')/(s'*s);
   plot(x(1),x(2),'ro','linewidth',2); axis square; drawnow; hold on; pause(.1)
   ns = norm(s); nf = norm(f); [ns nf];
   if ns < 1.e-12; break; end;

end;

[ x(1) x(2) f(1) f(2) ], xi=x(1);eta=x(2); 

p1 = xi/a2; q1 = eta/a2;

size_ratio = exp(pi/q1),      % Moffatt (3.11a)
strg_ratio = exp(pi*p1/q1),   % Moffatt (3.12)

