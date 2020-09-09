%% Purely Convection (Dirichlet Boundaty Condition)
close all;clear;clc;
%Physical Parameters
u = 0.5 ; a = 0.01; %a is alpha
%Parameters of RKW3
b1 = 1/4; b3 = 3/4; a21 = 8/15; a31 = 1/4; a32 = 5/12;
%Spatial Discretization
dx = 0.001;
N = 1/dx+1;
x = linspace(0,1,N);
%Time Discretization; For stable solution of purely convection equation, we
%need CFL = u(delta_t)/(delta_x) <= 1.73; Since delta_x is defined, this condtion can
%be Applied to determin maximum time step can take;
%delta_t <= %1.73*delta_x/u;
dt_max = 1.73*dx/u;
dt = 0.0001;
t0 = 0; tf = 2;
t = t0:dt:tf;
nof_t = 2/dt;
%Define initial condition;
T(:,1) = sin(pi*x).*(cos(5*pi*x)+sin(20*pi*x));%First Column; Initial Condition;
%Right hand side function; N-2 Equations;
%T_right from 3 to N, T_left from 1 to N-2;
RHS = @(T_right,T_left) -u*(T_right-T_left)/(2*dx);
k1 = zeros(N,1); k2 = zeros(N,1); k3 = zeros(N,1); %Preallocation; 
for i = 1:nof_t
    %Apply RKW3 to march interier spatial grid points from 2 to N-1 in time.
    %k1, k2, k3 are column vectors
    %Nothing happened for Boundary Points
    k1(1,1) = 0; k1(N,1) = 0;
    k1(2:N-1,1) = dt*RHS(T(3:N,i),T(1:N-2,i)); 
    
    k2(1,1) = 0; k2(N,1) = 0;
    k2(2:N-1,1) = dt*RHS(T(3:N,i)+a21*k1(3:N,1),T(1:N-2,i)+a21*k1(1:N-2,1)); 
    
    k3(1,1) = 0; k3(N,1) = 0;
    k3(2:N-1) = dt*RHS(T(3:N,i)+a31*k1(3:N,1)+a32*k2(3:N,1),T(1:N-2,i)+a31*k1(1:N-2,1)+a32*k2(1:N-2,1));
    
    T(:,i+1) = T(:,i) + b1*k1 + b3*k3;
end
%plot(x,T(:,1),x,T(:,1+0.5/dt),'-.',x,T(:,1+1/dt),':',x,T(:,end),'--');
figure
subplot(2,2,1); plot(x,T(:,1));title('t = 0s');
subplot(2,2,2); plot(x,T(:,1+0.5/dt)); title('t = 0.5s');
subplot(2,2,3); plot(x,T(:,1+1/dt)); title('t = 1s');
subplot(2,2,4); plot(x,T(:,end)); title('t = 2s');
%% Purely Convection (First Order Backward Difference for Right End Point)
close all;clear;clc;
%Physical Parameters
u = 0.5 ; a = 0.01; %a is alpha
%Parameters of RKW3
b1 = 1/4; b3 = 3/4; a21 = 8/15; a31 = 1/4; a32 = 5/12;
%Spatial Discretization
dx = 0.001;
N = 1/dx+1;
x = linspace(0,1,N);
dt = 0.0001;
t0 = 0; tf = 2;
t = t0:dt:tf;
nof_t = 2/dt;
T2(:,1) = sin(pi*x).*(cos(5*pi*x)+sin(20*pi*x));%First Column; Initial Condition;
%%All of this above same as previous problem;

%Right hand side function; N-2 Equations;
%T_right from 3 to N, T_left from 1 to N-2;
RHS = @(T_right,T_left) -u*(T_right-T_left)/(2*dx); %Central Differene for Interior Points
RHS2 = @(T_middle,T_left) -u*(T_middle-T_left)/dx; %Backward Differenc for Right Boundary Point
for i = 1:nof_t
    %Apply RKW3 to march interior spatial grid points from 2 to N;
    %Point 2 to N-1 by central difference; Point N by backward difference.
    
    %k1, k2, k3 are column vectors,returns increment in T;
    %Nothing happend for left boundary point
    k1(1,1) = 0; k1(N,1) = dt*RHS2(T2(N,i),T2(N-1,i)); %Boundary Points
    k1(2:N-1,1) = dt*RHS(T2(3:N,i),T2(1:N-2,i)); %Interior Points
    
    k2(1,1) = 0; k2(N,1) = dt*RHS2(T2(N,i)+k1(N,1),T2(N-1,i)+k1(N-1,1)); %Boundary Points
    k2(2:N-1,1) = dt*RHS(T2(3:N,i)+a21*k1(3:N,1),T2(1:N-2,i)+a21*k1(1:N-2,1));
    
    k3(1,1) = 0; k3(N,1) = dt*RHS2(T2(N,i)+k1(N,1)+k2(N,1),T2(N-1,i)+k1(N-1,1)+k2(N-1,1)); %Boundary Points
    k3(2:N-1,1) = dt*RHS(T2(3:N,i)+a31*k1(3:N,1)+a32*k2(3:N,1),T2(1:N-2,i)+a31*k1(1:N-2,1)+a32*k2(1:N-2,1));
    
    T2(:,i+1) = T2(:,i) + b1*k1 + b3*k3;
end
%plot(x,T(:,1),x,T(:,1+0.5/dt),'-.',x,T(:,1+1/dt),':',x,T(:,end),'--');
figure
subplot(2,2,1); plot(x,T2(:,1));title('t = 0s');
subplot(2,2,2); plot(x,T2(:,1+0.5/dt)); title('t = 0.5s');
subplot(2,2,3); plot(x,T2(:,1+1/dt)); title('t = 1s');
subplot(2,2,4); plot(x,T2(:,end)); title('t = 2s');
%% Purely Diffusion
close all;clear;clc
%Physical Parameters
u = 0.5 ; a = 0.01; %a is alpha
%Parameters of RKW3bv
dx = 0.001;
b1 = 1/4; b3 = 3/4; a21 = 8/15; a31 = 1/4; a32 = 5/12;
dt_max = (0.63*dx^2)/a; %D = a*dt/dx^2
dt = 1e-5;
N = 1/dx+1;
x = linspace(0,1,N); %Spatial Discretization
t0 = 0 ; tf = 10;
nof_t = round(tf/dt);
t = t0:dt:tf; %Time Discretization
%Initial Condition
T(:,1) = sin(pi*x).*(cos(5*pi*x)+cos(20*pi*x)); T(N,1) = 1;
%RHS function of dT/dt;
RHS = @(T_right,T_middle,T_left) a*(T_right-2*T_middle+T_left)/dx^2;
k1 = zeros(N,1); k2 = zeros(N,1); k3 = zeros(N,1); %Preallocation; 
for i = 1:nof_t %Marching Time Forward
    k1(1,1) = 0 ; k1(N,1) = 0; %Two Boundary Point;
    k1(2:N-1,1) = dt*RHS(T(3:N,i),T(2:N-1,i),T(1:N-2,i));
    
    k2(1,1)=0 ; k2(N,1) = 0;
    k2(2:N-1,1) = dt*RHS(T(3:N,i)+a21*k1(3:N,1),T(2:N-1,i)+a21*k1(2:N-1,1),T(1:N-2,i)+a21*k1(1:N-2,1));
    
    k3(1,1) = 0 ; k3(N,1) = 0;
    k3(2:N-1,1) = dt*RHS(T(3:N,i)+a31*k2(3:N,1)+a32*k2(3:N,1),T(2:N-1,i)+a31*k2(2:N-1,1)+a32*k2(2:N-1,1),...
        T(1:N-2,i)+a31*k2(1:N-2,1)+a32*k2(1:N-2,1));
    
    T(:,i+1) = T(:,i) + b1*k1 + b3*k3;
end
%Store temperature at required time
t_plot = [1 1+round(0.025/dt) 1+round(0.1/dt) 1+round(0.5/dt) 1+round(1/dt) 1+round(10/dt)];
for j = 1:6
    T_plot(:,j) = T(:,t_plot(j));
end
figure
subplot(3,2,1); plot(x,T_plot(:,1)); title('t = 0s');
subplot(3,2,2); plot(x,T_plot(:,2)); title('t = 0.025s');
subplot(3,2,3); plot(x,T_plot(:,3)); title('t = 0.1s');
subplot(3,2,4); plot(x,T_plot(:,4)); title('t = 0.5s');
subplot(3,2,5); plot(x,T_plot(:,5)); title('t = 1s');
ylim([-3 3]);
subplot(3,2,6); plot(x,T_plot(:,6)); title('t = 10s');
save('Diffusion','x','T_plot')