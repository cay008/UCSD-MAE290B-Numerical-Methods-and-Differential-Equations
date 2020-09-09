%% 2(b) Solution by ode23s
close all;clear;clc
tspan = [0 1000];
C0 = [0.5 7.5e-5 0 0];%Initial Conditions
[t,C] = ode23s(@odefnc, tspan,C0);
plot(t,C(:,1),t,C(:,2),'-.',t,C(:,3),':',t,C(:,4),'--','LineWidth',1.5)
title('Concentration of Chemicals Over Time');
legend('Chemical B','Chemical A','Chemical D','Chemical P');
xlabel('Time t'); ylabel('Concentration');


%% 2(b) Solution by RK4
close all;clear;clc
tic
%Define system of equations, f1 f2 f3 and f4 are RHS the given system
%y1 is chemical B, y2 is chemical A, y3 is chemical D, y4 is chemical P
k1 = 2.1e3; k2 = 5e-3; k3 = 18;
f1 = @(y1,y2,y3) -k1*y1*y2+k2*y3; %Function of B,A,D
f2 = @(y1,y2,y3) -k1*y1*y2+(k2+k3)*y3; %Function of B,A,D
f3 = @(y1,y2,y3) k1*y1*y2-(k2+k3)*y3; %Function of B,A,D
f4 = @(y3) k3*y3; %Function of D only
h = 0.0001; %Time step
t = 0:h:500; %Discretize time space; Equilibrium will be established before 500s.
%Initial Conditions
C(1,1) = 0.5; C(1,2) = 7.5e-5 ; C(1,3) = 0 ; C(1,4) = 0;
%Coefficients of RK4
a21 = 1/2 ; a32 = 1/2 ; a43 = 1;
c2 = 1/2 ; c3 = 1/2 ; c4 = 1;
b1 = 1/6 ; b2 = 1/3 ; b3 = 1/3 ; b4 = 1/6;
%Implement of RK4 to solve this system of 1st order ODEs
for j = 1:size(t,2)-1
    %Subscript ij represent ith Chemical, jth substep in RK4
    kB1 = h*f1(C(j,1) , C(j,2) , C(j,3));
    kA1 = h*f2(C(j,1) , C(j,2) , C(j,3));
    kD1 = h*f3(C(j,1) , C(j,2) , C(j,3));
    kP1 = h*f4(C(j,3) );
    
    kB2 = h*f1(C(j,1)+a21*kB1 , C(j,2)+a21*kA1 , C(j,3)+a21*kD1);
    kA2 = h*f2(C(j,1)+a21*kB1 , C(j,2)+a21*kA1 , C(j,3)+a21*kD1);
    kD2 = h*f3(C(j,1)+a21*kB1 , C(j,2)+a21*kA1 , C(j,3)+a21*kD1);
    kP2 = h*f4(C(j,3)+a21*kD1 );
    
    kB3 = h*f1(C(j,1)+a32*kB2 , C(j,2)+a32*kA2 , C(j,3)+a32*kD2);
    kA3 = h*f2(C(j,1)+a32*kB2 , C(j,2)+a32*kA2 , C(j,3)+a32*kD2);
    kD3 = h*f3(C(j,1)+a32*kB2 , C(j,2)+a32*kA2 , C(j,3)+a32*kD2);
    kP3 = h*f4(C(j,3)+a32*kD2 );
     
    kB4 = h*f1(C(j,1)+a43*kB3 , C(j,2)+a43*kA3 , C(j,3)+a43*kD3);
    kA4 = h*f2(C(j,1)+a43*kB3 , C(j,2)+a43*kA3 , C(j,3)+a43*kD3);
    kD4 = h*f3(C(j,1)+a43*kB3 , C(j,2)+a43*kA3 , C(j,3)+a43*kD3);
    kP4 = h*f4(C(j,3)+a43*kD3 );
    
    C(j+1,1) = C(j,1) + b1*kB1 + b2*kB2 +b3*kB3 + b4*kB4 ;
    C(j+1,2) = C(j,2) + b1*kA1 + b2*kA2 +b3*kA3 + b4*kA4 ;
    C(j+1,3) = C(j,3) + b1*kD1 + b2*kD2 +b3*kD3 + b4*kD4 ;
    C(j+1,4) = C(j,4) + b1*kP1 + b2*kP2 +b3*kP3 + b4*kP4 ;
end
plot(t,C(:,1),t,C(:,2),'-.',t,C(:,3),':',t,C(:,4),'--','LineWidth',1.5)
title({'Chemical Concentration Over Time','RK4'});
xlabel('Time t');ylabel('Chemical Concentration');
legend('Chemical B','Chemical A','Chemical D','Chemical P');
toc


    
%% 3(a) Solution by Shooting Method
close all; clear;clc
%Equations
a = @(x) -(x+3)./(x+1);
b = @(x) (x+3)./(x+1).^2;
c = @(x) 2*(x+1)+3*b(x);
f = @(z) z ; %RHS of first ode; function of z only
g = @(x,y,z) -a(x).*z -b(x).*y +c(x) ; %RHS of second ode; function of x y z

%All essentials for solution
h = 0.01; %Time step
x = 0:h:2;


y1(1) = 5 ; %First boundary conditoon @ x = 0;
y2(1) = 5;

z1(1) = 1 ; %First Guess of BC
z2(1) = 2 ; %Second Guess of BC

%Coefficients of RK
a21 = 1/2 ; a32 = 1/2 ; a43 = 1;
c2 = 1/2 ; c3 = 1/2 ; c4 = 1;
b1 = 1/6 ; b2 = 1/3 ; b3 = 1/3 ; b4 = 1/6;

%Implement of RK4
for j = 1:size(x,2)-1
    %k for increment in y, subscript ij represent ith guess, jth substep in RK4
    %l for increment in z
    
    %RK4 Corresponding to first guess
    
    k11 = h*f(z1(j));
    l11 = h*g(x(j) , y1(j) , z1(j));
    
    k12 = h*( f(z1(j) + a21*l11) ) ;
    l12 = h*( g( x(j)+c2*h , y1(j)+a21*k11 , z1(j)+a21*l11 ) );
     
    k13 = h*( f(z1(j)+a32*l12) ) ;
    l13 = h*( g(x(j)+c3*h , y1(j)+a32*k12 , z1(j)+a32*l12) ) ;
     
    k14 = h*( f(z1(j)) + a43*l13) ;
    l14 = h*( g(x(j)+c4*h , y1(j)+a43*k13 , z1(j)+a43*l13) ); 
     
    y1(j+1) = y1(j) + b1*k11 + b2*k12 + b3*k13 + b4*k14;%First Solution
    z1(j+1) = z1(j) + b1*l11 + b2*l12 + b3*l13 + b4*l14;
    
    %RK4 Corresponding to second guess
    k21 = h*f(z2(j));
    l21 = h*g(x(j),y2(j),z2(j));
    
    k22 = h*( f(z2(j)+a21*l21) ) ;
    l22 = h*( g(x(j)+c2*h,y2(j)+a21*k21,z2(j)+a21*l21) );
     
    k23 = h*( f(z2(j)+a32*l22) ) ;
    l23 = h*( g(x(j)+c3*h,y2(j)+a32*k22,z2(j)+a32*l22) ) ;
     
    k24 = h*( f(z2(j)) + a43*l23) ;
    l24 = h*( g(x(j)+c4*h,y2(j)+a43*k23,z2(j)+a43*l23) ); 
     
    y2(j+1) = y2(j) + b1*k21 + b2*k22 + b3*k23 + b4*k24;%Second Solution
    z2(j+1) = z2(j) + b1*l21 + b2*l22 + b3*l23 + b4*l24;
    
    j=j+1
end
%Find the linear combinations of this two solutions to satisfy boundary
%condition T_1(x=L)=4, coeffij represents ith BC with jth coefficient
coeff11 = (4-y2(1,end))/(y1(1,end)-y2(1,end));
coeff12 = (y1(1,end)-4)/(y1(1,end)-y2(1,end));
y = coeff11*y1 + coeff12*y2;
figure; plot(x,y,'LineWidth',1.5); 
title({'Temperature Distribution Along the Bar','Dirichlet Boundary Condition'})
xlabel('x Position'); ylabel('Temperature');

coeff21 = (0-z2(1,end))/(z1(1,end)-z2(1,end));
coeff22 = (z1(1,end)-0)/(z1(1,end)-z2(1,end));
y_2ndBC = coeff21*y1 + coeff22*y2;
figure
plot(x,y_2ndBC,'LineWidth',1.5); 
title({'Temperature Distribution Along the Bar','Neumann Boundary Condition'})
xlabel('x Position'); ylabel('Temperature');



%% 3(b) Solution by Thomas Algorithm ##Need Run Last Section First##
%Keep my solutions from shooting methods, since I need to compare them with
%solution from direct methods.
clearvars -except x y;
x_shooting = x;
y_shooting = y;

x = linspace(0,2,21);%Discretize the space to 20 elements,total 21 Nodes
h = 2/(21-1); %Step Size
%Find all the coefficients in discretized equation
for i = 1:21 
    a(i) = -(x(i)+3)/((x(i)+1));
    b(i) = (x(i)+3)/(x(i)+1).^2;
    f(i) = 2.*(x(i)+1)+3*b(i);
    alpha(i) = 1-a(i)*h/2;
    beta(i) = b(i)*h^2-2;
    gama(i) = 1+a(i)*h/2;
end

%Construct the linear system to be solved, total 19 Equations; 
%First and last equation need to be treated seperately since they contain 
%known values from boundary conditions
A(1,1) = beta(2)      ; A(1,2) = gama(2);
A(19,18) = alpha(20)  ; A(19,19) = beta(20);

delta(1,1) = f(2)*h^2-5*alpha(2);
delta(19,1) = f(20)*h^2-4*gama(20);

%Then deal with remaining 17 Equations
for j = 2:18 
    A(j,j) = beta(j+1);
    A(j,j-1) = alpha(j+1);
    A(j,j+1) = gama(j+1);
    delta(j,1) = f(j+1)*h^2;
end
%Until now we finally obtained the linear system AT=delta (Ax=b)

s1 = A\delta

y_2ndBC = thomas(A,delta)
temperature(1) = 4; temperature(2:20)=y_2ndBC; temperature(21)=4
temp2(1) = 5; temp2(2:20)=s1 ; temp2(21)=4
figure
plot(x,temp2,'*','LineWidth',1.5);hold on; 
plot(x_shooting,y_shooting);hold off
title('Comparison Between Shooting Method and Direct Method')
xlabel('x Position');ylabel('Temperature');
legend('Direct Method','Shooting Method');

%% System of equations which need to be solved in 2(b)
function dCdt = odefnc(t,C)
k1 = 2.1e3; k2 = 5e-3; k3 = 18
dCdt = [-k1*C(1)*C(2)+k2*C(3);-k1*C(1)*C(2)+(k2+k3)*C(3);k1*C(2)*C(1)-...
    (k2+k3)*C(3); k3*C(3)]
%C(1) = C_B ; C(2) = C_A ; C(3) = C_D ; C(4) = C_P
end


%% Thomas Algorithm for 3((b)
function solution = thomas(A,b)
s = size(A,1); %Number of rows
%Forward Elimination of subdiagonal elements;
for j = 1:s-1;		

    M(j) = -A(j+1,j)/A(j,j);%Multiplier
    A(j+1,j) = A(j,j)*M(j) + A(j+1,j);
    A(j+1,j+1) = A(j+1,j+1)+A(j,j+1)*M(j);%New next diagonal term
    b(j+1) = b(j+1)+M(j)*b(j);
    %No need to modify subdiagonal terms since they are not used in the
    %following steps
    j=j+1
end
%Backward substitution for solution;
solution(s) = b(s)/A(s,s); %The Last element of solution
for k = s-1:-1:1
    solution(k) = (b(k)-A(k,k+1)*solution(k+1))/A(k,k);%Update b vector with solutions  gradually
end
end
