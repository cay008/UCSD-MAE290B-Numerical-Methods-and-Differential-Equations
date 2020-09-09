%% Problem1
close all;clear;clc
tic
%Discretization of Equation
dx=0.01;x0=0;xf=20; x=x0:dx:xf; N=size(x,2);
dy = 0.0005; 
y0=0;yf=35; y=y0:dy:yf; noft=size(y,2); %y is time like variable
%Parameters of RKW3
b1 = 1/4; b3 = 3/4; a21 = 8/15; a31 = 1/4; a32 = 5/12;
%Second Order Central Scheme
RHS = @(phi_right,phi_middle,phi_left) (-i/20)*(phi_right-2*phi_middle+phi_left)/dx^2;
%Initializations
phi = zeros(N,noft);
phi_inter = zeros(N,1);
k = zeros(N,1);
%Initial Condition
phi(:,1) = exp(-(x-5).^2/4)+exp(-(x-15).^2/4+10*i*x);
%Apply RKW3 to march phi forward in y
for i = 1:noft
    phi(1,i) = 0; phi(end,i) = 0; 
    phi(:,i+1) = phi(:,i); %phi_n+1 = phi_n so far;
    %k1
    k(1,1)=0; k(end,1)= 0;
    k(2:N-1,1) = dy*RHS(phi(3:N,i),phi(2:N-1,i),phi(1:N-2,i));
    %March solution forward for one substep
    phi(:,i+1) = phi(:,i+1) + b1*k; %phi_n+1 = phi_n + b1*k1 so far;
    %Increment value of u that needed for k2
    phi_inter = phi(:,i) + a21*k; 
    %k2
    k(2:N-1,1) = dy*RHS(phi_inter(3:N,1),phi_inter(2:N-1,1),phi_inter(1:N-2,1));
    %Next value to be plugged into RHS function
    phi_inter = phi(:,i+1) + a32*k; %phi_n + b1(same as a31)*k1 + a32*k2
    %k3
    k(2:N-1,1) = dy*RHS(phi_inter(3:N,1),phi_inter(2:N-1,1),phi_inter(1:N-2,1));
    %Solution of next time step;
    phi(:,i+1) = phi(:,i+1) + b3*k;
end
%Store Part of Result to a new matrix for contour plot
P = phi(:,1:100:noft);
P = (abs(P)).^2;
Plot = pcolor(P);shading flat;colorbar;
set(gca,'XTickLabel',[5 10 15 20 25 30 35]);
set(gca,'YTickLabel',linspace(2,20,10));
xlabel('Time Like Variable : y');ylabel('x');
title('Contour Plot')
toc
saveas(gcf,'Laser Beam','bmp');
%% Problem2 (a) Second Order Central Scheme
close all;clear;clc;
%Parameters of Physical Problem
a = 1/2;
%Parameters of RKW3
b1 = 1/4; b3 = 3/4; a21 = 8/15; a31 = 1/4; a32 = 5/12;
%Discretization
dx=0.01; x0=0; xf=5;x=x0:dx:xf; N=size(x,2);
xc = 1 + 2/dx; %Critical Point;
dt=0.001; t0=0; tf=10; t=t0:dt:tf; nof_t=size(t,2);
%Initializations
u = zeros(N,nof_t);
u_inter = zeros(N,1);
k = zeros(N,1);
%RHS Function of Second order Central Scheme for interior points
RHS_C = @(T_right,T_left) -(1/2)*(T_right.^2-T_left.^2)/(2*dx); 
%RHS Function of First order Backward Scheme for interior points
RHS_B = @(T_middle,T_left) -(1/2)*(T_middle.^2-T_left.^2)/dx;
%Initial Condition
u(1:xc,1) = exp(-((x(1:xc))-1).^2/0.18);
u(xc:end,1) = zeros(size(u(xc:end,1),2),1);
for i = 1:nof_t
    %Visualize how this wave will propagate
    %plot(x,u(:,i));ylim([-0.5 1.5]);pause(0.01);
    u(:,i+1) = u(:,i); %u_n+1 = u_n so far;
    %k1
    k(1,1)=0; k(end,1)= dt*RHS_B(u(N,i),u(N-1,i));
    k(2:N-1,1) = dt*RHS_C(u(3:N,i),u(1:N-2,i));
    %March solution forward for one substep
    u(:,i+1) = u(:,i+1) + b1*k; %u_n+1 = u_n + b1*k1 so far;
    %Increment value of u that needed for k2
    u_inter = u(:,i) + a21*k; 
    %k2
    k(2:N-1,1) = dt*RHS_C(u_inter(3:N,1),u_inter(1:N-2,1));
    k(N,1) = dt*RHS_B(u_inter(N,1),u_inter(N-1,1));
    %Next value to be plugged into RHS function
    u_inter = u(:,i+1) + a32*k; %u_n + b1(same as a31)*k1 + a32*k2
    %k3
    k(2:N-1,1) = dt*RHS_C(u_inter(3:N,1),u_inter(1:N-2,1));
    %Solution of next time step;
    u(:,i+1) = u(:,i+1) + b3*k;
end
Plot = figure; 
sb1 = subplot(3,2,1); plot(x,u(:,1)); sb2 = subplot(3,2,2); plot(x,u(:,1+0.5/dt));
sb3 = subplot(3,2,3); plot(x,u(:,1+1/dt)); sb4 = subplot(3,2,4); plot(x,u(:,1+1.5/dt));
sb5 = subplot(3,2,5); plot(x,u(:,1+5/dt)); sb6 = subplot(3,2,6); plot(x,u(:,end));
ylim([-0.5,1.5]);
linkaxes([sb1,sb2,sb3,sb4,sb5,sb6],'y');
x_all = get(findobj(Plot,'Type','Axes'),'Xlabel');
y_all = get(findobj(Plot,'Type','Axes'),'YLabel');
title_all = get(findobj(Plot,'Type','Axes'),'Title')
set([x_all{:}],'String','x');
set([y_all{:}],'String','u');
set([title_all{:}],'String','Evolution of Wave');
%% Problem2 (b) First Order Upwind Scheme
close all;clear;clc
%Parameters of Physical Problem
a = 1/2;
%Parameters of RKW3
b1 = 1/4; b3 = 3/4; a21 = 8/15; a31 = 1/4; a32 = 5/12;
%Discretization
dx=0.01; x0=0; xf=5;x=x0:dx:xf; N=size(x,2);
xc = 1 + 2/dx; %Critical Point;
dt=0.001; t0=0; tf=10; t=t0:dt:tf; nof_t=size(t,2);
%Initializations
u = zeros(N,nof_t);
u_inter = zeros(N,1);
k = zeros(N,1);
%RHS Function of Upwind Scheme
RHS_U = @(T_middle,T_left) -(1/2)*(T_middle.^2-T_left.^2)/dx
%Initial Condition
u(1:xc,1) = exp(-((x(1:xc))-1).^2/0.18);
u(xc:end,1) = zeros(size(u(xc:end,1),2),1);
for i = 1:nof_t
    %Visualize how this wave will propagate
    %plot(x,u(:,i));ylim([-0.5 1.5]);pause(0.01);
    u(:,i+1) = u(:,i); %u_n+1 = u_n so far;
    %k1
    k(1,1)=0; k(2:N,1) = dt*RHS_U(u(2:N,i),u(1:N-1,i));
    %March solution forward for one substep
    u(:,i+1) = u(:,i+1) + b1*k; %u_n+1 = u_n + b1*k1 so far;
    %Increment value of u that needed for k2
    u_inter = u(:,i) + a21*k; 
    %k2
    k(2:N,1) = dt*RHS_U(u_inter(2:N,1),u_inter(1:N-1,1));
    %Next value to be plugged into RHS function
    u_inter = u(:,i+1) + a32*k; %u_n + b1(same as a31)*k1 + a32*k2
    %k3
    k(2:N,1) = dt*RHS_U(u_inter(2:N,1),u_inter(1:N-1,1));
    %Solution of next time step;
    u(:,i+1) = u(:,i+1) + b3*k;
end
Plot = figure; 
sb1 = subplot(3,2,1); plot(x,u(:,1)); sb2 = subplot(3,2,2); plot(x,u(:,1+0.5/dt));
sb3 = subplot(3,2,3); plot(x,u(:,1+1/dt)); sb4 = subplot(3,2,4); plot(x,u(:,1+1.5/dt));
sb5 = subplot(3,2,5); plot(x,u(:,1+5/dt)); sb6 = subplot(3,2,6); plot(x,u(:,end));
ylim([-0.5,1.5]);
linkaxes([sb1,sb2,sb3,sb4,sb5,sb6],'y');
x_all = get(findobj(Plot,'Type','Axes'),'Xlabel');
y_all = get(findobj(Plot,'Type','Axes'),'YLabel');
title_all = get(findobj(Plot,'Type','Axes'),'Title')
set([x_all{:}],'String','x');
set([y_all{:}],'String','u');
set([title_all{:}],'String','Evolution of Wave');