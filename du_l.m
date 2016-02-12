function [ d ] = du_l( x,u,params )
%derivative of loss for cart pole problem wrt to u evaluated
%at state vector x, control vector u
%params = structure of parameters

d(1) = 2*params.Fcost*u(1);
d(2) = 2*params.Tcost*u(2);
d=d';

end

