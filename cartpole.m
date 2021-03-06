clear all;
%main script for doing gradient-based optimal control on the
%underactuated cart-pole system. run me!

params.dt = 0.005; %euler integration step length, in seconds
params.m1 = 0.5; %mass of cart, in kg
params.m2 = 20.5; %mass at the end of the pole, in kg
params.l = 1; %length of pole, in meters
params.g = 9.81; %gravity, in m/s^2
params.mu = 1; %viscous friction coefficient
params.nsteps = 400; %total number of timesteps
params.T = params.nsteps-1; %number of equality constraints/number of transitions
params.Tcost = 50; %cost weight on the cart-pole joint torque
params.Fcost = 1e-4; %cost weight on F, the force applied to the cart
params.xcost = 100; %cost weight for not being at theta=-3.14 rads
muscale = 0.9; %amount to scale mu on each outer iteration/warm restart
Tcostscale = 4; %amount to scale the cost weight on the cart-pole joint torque for each outer iteration
x0 = [0;0;0;0]; %initial state

best_us = (rand((params.nsteps-1)*2,1)-0.5)*50; %initial open loop controls
outeriter = 5; %number of outer iterations/warm restarts

for i = 1:outeriter
    params.Tcost = params.Tcost*Tcostscale;
    params.mu = params.mu*muscale;
    fprintf('On outer iteration %d/%d\n',i,outeriter);
    fprintf('Mu (viscous damping coefficient) and weight on T control are:\n')
    fprintf('Mu: %f \nTcost: %f \n',params.mu,params.Tcost);
    
    fun = @(x) cartpole_grad(x,params,x0);
    options = optimoptions('fminunc','MaxIter',400,'GradObj','on',...
        'Algorithm','quasi-newton','PlotFcns',@optimplotfval);
    [best_us,fval] = fminunc(fun,best_us,options);
    shaped_us = reshape(best_us,[params.nsteps-1,2])'; %put back to 2xT vector
    fprintf('Norm of cart-pole joint torque is (this should go towards zero): %f\n',norm(shaped_us(2,:)))
end
%%
%Turn off the cart-pole joint torque and plot the result
shaped_us(2,:) = 0;
[loss,xs] = sim_loss(x0,shaped_us,params);
times = (1:1:params.nsteps)*params.dt;
figure;
plot(times,xs(2,:));
xlabel('Time (s)');
ylabel('Theta (rad)');