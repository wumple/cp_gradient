function [ loss,grad ] = cartpole_grad(us,params,x0)
    %computes the gradient of the lagrangian with respect to the controls
    %'us' should be a m x T matrix, m = dimension of control vector, T =
    %number of controls
    %x0 is the starting vector
    %params is a structure of parameters, see main script cartpole.m
    
    us = reshape(us,[params.nsteps-1,2])'; %for compatibility with Matlab's optimizer...
    
    %first simulate the system forward to get the loss and the state
    %vectors
    loss = 0;
    xs = zeros(4,params.nsteps);
    xs(:,1) = x0;
    nsteps = params.nsteps;
    T = params.T;
    for i=1:nsteps-1
        loss = loss + loss_cp(xs(:,i),us(:,i),params);
        xs(:,i+1) = step_cp(xs(:,i),us(:,i),params);
    end

    loss = loss + loss_cp(xs(:,nsteps),zeros(2,1),params);
    
    %only run this if a gradient is requested:
    if(nargout > 1)
    
        %calculate the lambdas with the recurrence relation...
    
        lambdas = zeros(4,T);
        lambdas(:,T) = -dx_l(xs(:,end),0,params);
        counter = T-1;

        while counter > 0
            lambdas(:,counter) = -(dx_l(xs(:,counter+1),us(:,counter+1),params)-...
                jac_f(xs(:,counter+1),us(:,counter+1),params)*lambdas(:,counter+1));
            counter = counter-1;
        end
    
        %now use the lambdas to calculate the actual gradient
        us_grad = zeros(size(us));    
        for i=1:T
            dul = du_l(xs(:,i),us(:,i),params);
            duf = du_f(xs(:,i),us(:,i),params);
            for p = 1:size(us,1)
                us_grad(p,i) = dul(p) - sum(duf(p,:)'.*lambdas(:,i));
            end
        end
    
        %for optimizer compatibility...
        grad = reshape(us_grad',[2*(params.nsteps-1),1]);

    end

