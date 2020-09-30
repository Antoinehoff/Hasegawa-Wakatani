function [ time, Napj ] = RK4(f, time, f0)
%RK4 explicit ODE solver
Napj      = zeros(numel(f0),numel(time));
Napj(:,1) = f0;
for it = 1:(numel(time)-1)
    dt  = time(it+1) - time(it);
    Y   = Napj(:,it);
    k_1 = f(Y);
    k_2 = f(Y + 0.5*dt*k_1);
    k_3 = f(Y + 0.5*dt*k_2);
    k_4 = f(Y +     dt*k_3);    
    
    Napj(:,it+1) = Y + dt*(1./6.)*(k_1 + 2.*k_2 + 2.*k_3 + k_4);
end

Napj = transpose(Napj);
time = transpose(time);
end

