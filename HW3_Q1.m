%% 1(c)
close all; clear; clc;
dt = 0.0001;
dx = 1/40;
[x,T_pt0001]=FE(0,100,41,dt);
figure
plot(x,T_pt0001(:,1+1/dt),x,T_pt0001(:,1+10/dt),'--',x,T_pt0001(:,end),'-.');
legend('T=1s','T=10s','T=100s');
xlabel('Spatial Domain'); ylabel('Temperature');title('Temperature Distribution in the domain');
grid on
%% 1(d)
close all;clear;clc
dx = 0.001; dt = 0.0001 ; N = 1+1/dx;
[x,T] = FE(0,1,N,dt);
T_abs = T(1+0.55/dx,end);%This the standard temperature for error calculation;
N = 41:40:1001
delta_x = 1./(N-1);%Changing delta_x
tic
for i = 1:max(size(N))
    [x,T] = FE(0,1,N(i),dt); 
    error(i) = abs(T(1+round(0.55/delta_x(i)),end)-T_abs);
end
toc
figure
plot(delta_x,error);
xlabel('\Deltax');ylabel('Error');title('Effect of Changing \Deltax');
%% 1(e)
close all;clear;clc
dx = 0.001; dt = 0.0001 ; N = 1+1/dx;
[x,T] = FE(0,1,N,dt);
T_abs = T(1+0.55/dx,end);%This the standard temperature for error calculation;
delta_t = linspace(0.01,0.0001,20);%Changing delta t;
tic
for i = 1:max(size(delta_t))
    [x,T] = FE(0,1,N,delta_t(i)); 
    error(i) = abs(T(1+round(0.55/dx),end)-T_abs);
end
toc
plot(delta_t,error);
xlabel('\Deltat');ylabel('Error');title('Effect of Changing \Deltat');
%% 2(d) (i) Dirichlet Boundary Conditon
close all;clear;clc;
%Physical Parameters
u = 0.5 ; alpha = 0.01;
%Parameters of RKW3
b1 = 1/4; b3 = 3/4; a21 = 8/15; a31 = 1/4; a32 = 5/12;
%Spatial Discretization
dx = 0.01;
N = 1/0.01+1;
x = linspace(0,1,N);
%Time Discretization; For stable solution of purely convection equation, we
%need CFL = u(delta_t)/(delta_x) <= 1.73; Since delta_x is defined, this condtion can
%be Applied to determin maximum time step can take; 
%delta_t <= %1.73*delta_x/u;
dt_max = 1.73*dx/u; %dt_max = 0.0346
dt = 0.001;
t0 = 0; tf = 2;
t = t0:dt:tf;
nof_t = 2/dt;
%Define initial condition;
T(:,1) = sin(pi*x').*(cos(5*pi*x')+sin(20*pi*x'));
for i = 1:nof_t
    T(1,i+1) = 0; %Left boundary condition
    T(N,i+1) = 0; %Right Boundary Condition; Dirichlet;
    for j = 2:N-1 %Interier Nodes
        k1 = -(u*dt/dx^2)*(T(j+1)-2*T(j)+T(j-1));
        k3 = -(u*dt/dx^2)*(1+a31+a32)*(T(j+1)-2*T(j)+T(j-1));
        %March time forward for solution;
        T(j,i+1) = T(j,i+1) + b1*k1 + b3*k3;
    end
end
plot(x,T(:,1),x,T(:,1+0.5/dt),'-.',x,T(:,1+1/dt),':',x,T(:,end),'--')
title({'Solution of convection equaiton','Dirichlet Boundary Condition'});
xlabel('x');ylabel('T'); legend('t = 0','t = 0.5','t = 1','t = 2')
