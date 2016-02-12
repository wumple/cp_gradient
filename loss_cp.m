function [ l ] = loss_cp(x,u,params)

%cart pole loss for driving system to unstable equilibrium at theta (x(2))
% = -pi
%x = state vector
%u = control vector
%params = structure of parameters

l = params.xcost*(x(2)-pi)*(x(2)-pi) + params.Tcost*u(2)*u(2) + params.Fcost*u(1)*u(1);

end

