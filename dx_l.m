function [ d ] = dx_l( x,u,params )
%derivative of cart pole loss
%with respect to state vector x, control vector u
%params = structure of parameters

d = zeros(1,4);
d(1) = 0;
d(2) = 2*params.xcost*(x(2)-pi);
d(3) = 0;
d(4) = 0;
d = d';


end

