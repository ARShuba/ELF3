clear all
%%
tspan = 0:1:2000;
% initial condition
x_start1 = [33554.833280 56.5 0 0 0 0 0 0 0];
%calculating steady state for given initial condition 
% [t,x_time] = ode23s('dynamics',tspan,x_start);
[t1,x_time1] = ode23s('dynamics',tspan,x_start1);

% plot(t,x_time(:,3))
xlim([0,1600])
hold on
plot(t1,(x_time1(:,2)),'Linewidth',2)