function [ d ] = du_f( x,u,params )
%derivative of transition function f wrt controls u at state x, controls u
%params = structure of parameters

dt = params.dt;
m1 = params.m1;
m2 = params.m2;
M = m1+m2;
l = params.l;
g = params.g;
mu = params.mu; %viscous friction

d = zeros(2,4);

%derivs of first output wrt u1,u2
d(1,1) = 0;
d(2,1) = 0;

%derivs of second output wrt u1,u2
d(1,2) = 0;
d(2,2) = 0;

%derivs of third output wrt u1,u2
dx4dotdF = (1/(m2*l*cos(x(2))-M*(l+mu)/cos(x(2))));
d(1,3) = dt*(-(l+mu)/cos(x(2)))*dx4dotdF;
dx4dotdT = (-M/cos(x(2)))/(m2*l*cos(x(2))-M*(l+mu)/cos(x(2))); %T =u(2)
d(2,3) = dt*((1/(cos(x(2))))-(l+mu)*dx4dotdT/cos(x(2)));

%derivs of fourth output wrt u1,u2
d(1,4) = dt*dx4dotdF;
d(2,4) = dt*dx4dotdT; %d(4,1)





end

