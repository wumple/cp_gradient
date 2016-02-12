function [ loss,xs] = sim_loss( x0,us,params)
%simulates cart-pole system forward starting from x0
%with open loop controls 'us'
%params = structure of parameters


    loss = 0;
    xs = zeros(4,params.nsteps);
    xs(:,1) = x0;
    for i=1:params.nsteps-1
        loss = loss + loss_cp(xs(:,i),us(:,i),params);
        xs(:,i+1) = step_cp(xs(:,i),us(:,i),params);
    end

    loss = loss + loss_cp(xs(:,params.nsteps),zeros(2,1),params);

end

