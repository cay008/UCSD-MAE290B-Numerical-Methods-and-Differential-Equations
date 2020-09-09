%% (b) Matlab Code
close all;clear;clc
tic
%Physical Paramters
O = 30; alpha = 0.1; a = 3;
Q = @(x,y,t) 2.5*sin(3*pi.*x).*sin(4*pi.*y).*(1-exp(-a*t).*cos(O*t).*cos(O*t));

%Space Discretization
M = 81; x = linspace(0,1,81); dx = 1/(M-1);
N = 81; y = linspace(0,1,81); dy = 1/(N-1);
[X,Y] = meshgrid(x,y);

%Time
dt = 0.01;


%Operators
TD = gallery('tridiag',M-2,-1,2,-1);
I = speye(M-2,N-2); %Only 79 interior nodes need to be solved;
beta = alpha*dt/(2*dx^2);
gamma = alpha*dt/(2*dy^2);
A1 = I + beta*TD;   B1 = I - gamma*TD;
A2 = I + gamma*TD;  B2 = I - beta*TD;
E1 = I + gamma*TD;  F1 = I - beta*TD;
E2 = I + beta*TD;   F2 = 1 - gamma*TD;

%Initial Condition
T(:,:,1) = Q(X,Y,0); %Source Term;
T(:,[1,N],1) = 0; T([1,M],:,1) = 0; %Dirichlet Boundary Condition;

%Initialize the vector to Store Temperature at x = 0.55, y = 0.45; At this
%stage, T_plot start with two elements, and it will grow continuously
%during the loop of ADI;
X_plot = 1 + 0.55/dx; Y_plot = 1 + 0.45/dy;
T_plot = [T(X_plot,Y_plot),1];

%Implementation of ADI
i = 0;
t(1) = 0;
T_inter = zeros(M-2,N-2); %Temperature at half intermediate time step
T_store = zeros(M,N); %Store the obtained result of ADI
while 1
    i = i + 1;
    Q1 = (dt/2)*Q(X,Y,(i-1/2)*dt); %Source term at next half time step
    Q2 = (dt/2)*Q(X,Y,i*dt);       %Source term at next time step
    s = 1;                  %Parameter used to execute two scheme turn by turn
    
    if s == 1
        b = B1*T(2:M-1,2:N-1,i) + Q1(2:M-1,2:N-1);
        %Since A1 act on the rows of T (x direction),  then the returned
        %soltion of this for loop is the transpose matrix of the actual one
        T_inter = thomas(A1,b);
        %B2 act on the rows of actual T_n+1/2, same effect as acting on
        %columns of T_inter from the last loop.
        b = B2*T_inter + Q2(2:M-1,2:N-1);
        T_inter = thomas(A2,b);
        %Return T_n+1
        T_store(2:M-1,2:N-1) = T_inter;
        s = s - 1;
        
    else %Reverse the scheme of x and y
        %F1 should act on the rows of T_n, the effect of which is same as
        %F1 acting on the columns of T_n'
        b = F1*T(2:M-1,2:N-1,i).' + Q1(2:M-1,2:N-1);
        T_inter = thomas(E1,b);
        %Return T_n+1/2
        b = F2*T_inter + Q2(floor(2:M-1.2:N-1));
        %E2 should act on the rows of T_n+1, the effect of which is same as
        %E2 acting on the columns of (T_n+1/2)'
        T_inter = thomas(E2,b);
        %Return (T_n+1)'
        T_inter = T_inter.';
        T_store(2:M-1,2:N-1) = T_inter;
        s = s + 1;
    end
    %Update Solution
    T(:,:,i+1) = T_store;
    %Solution at the next time step at x = 0.55, y = 0.45
    T_plot(i+1) = T(X_plot,Y_plot,i+1);
    t(i+1) = dt*i; %Time vector
    %Convergence Criteria
    if abs(T_plot(i+1)-T_plot(i)) <= eps
        break
    end
end
%Visualize steady state solution
figure; contourf(T(:,:,end));colorbar;
set(gca,'XTickLabel',[dx*10 dx*20 dx*30 dx*40 dx*50 dx*60 dx*70 dx*80]);
set(gca,'YTickLabel',[dx*10 dx*20 dx*30 dx*40 dx*50 dx*60 dx*70 dx*80]);
xlabel('x'); ylabel('y'); title('Contour of Steady State Solution');

figure; surf(T(:,:,end));colorbar;
set(gca,'XTickLabel',[0 0.5 1]);
set(gca,'YTickLabel',[0 0.5 1]);
xlabel('x'); ylabel('y'); zlabel('Temperature');
title('Contour of Steady State Solution');

%Time revolution of temperature at x = 0.55, y = 0.45;
figure; plot(t,T_plot)
xlabel('Time'); ylabel('Temperature');
title('Temperature Revolution');
%Zoom in plot, t range from t = 0 to t = 2;
figure; plot(t(1:2/dt+1),T_plot(1:2/dt+1));
xlabel('Time'); ylabel('Temperature');
title('Zoom in Plot of Detailed Temperature Revolution');
toc


%Try to find how source term looks like over time;
%tau = 0:0.001:1;
%S = zeros(81,81,size(tau,2));
%for z = 1:size(tau,2)
%S(:,:,z) = Q(X,Y,tau(z));
%S_plot(z) = S(X_plot,Y_plot,z);
%end
%figure;plot(tau,S_plot)


%% (d) Matlab Code
close all;clear;clc
%Since in this part, we want to figure out how Temperature at x = 0.55 and
%y = 0.45 will revolute with varying Omega and thermal diffusitivity alpha;
%we can loose the tolerance
tic
%Physical Paramters
O = linspace(1,30,10)
alpha = linspace(0.1,5,10)
a = 3;

%Space Discretization
M = 81; x = linspace(0,1,81); dx = 1/(M-1);
N = 81; y = linspace(0,1,81); dy = 1/(N-1);
[X,Y] = meshgrid(x,y);

%Time Step
dt = 0.001;
t0 = 0 ; tf = 2; t = t0:dt:tf;
nof_t = size(t,2);

%
omega_count = 0 ;
alpha_count = 0 ;
T_plot = zeros(size(O,2)*size(alpha,2),nof_t);

for ac = 1 : size(alpha,2)
    alpha_count = alpha_count + 1;
    figure;
    for oc = 1 : size(O,2)                %Loop over all Omega
        omega_count = omega_count + 1;
        
        Q = @(x,y,t) 2.5*sin(3*pi.*x).*sin(4*pi.*y).*(1-exp(-a*t).*cos(O(oc)*t).*cos(O(oc)*t));
        
        %Operators
        TD = gallery('tridiag',M-2,-1,2,-1);
        I = speye(M-2,N-2); %Only 79 interior nodes need to be solved;
        beta = alpha(ac)*dt/(2*dx^2);
        gamma = alpha(ac)*dt/(2*dy^2);
        A1 = I + beta*TD;   B1 = I - gamma*TD;
        A2 = I + gamma*TD;  B2 = I - beta*TD;
        E1 = I + gamma*TD;  F1 = I - beta*TD;
        E2 = I + beta*TD;   F2 = 1 - gamma*TD;
        
        %Initial Condition
        T(:,:,1) = Q(X,Y,0); %Source Term;
        T(:,[1,N],1) = 0; T([1,M],:,1) = 0; %Dirichlet Boundary Condition;
        
        %Initialze the matrix to Store all temperature revolution at x = 0.55
        %and y = 0.45.
        X_plot = 1 + 0.55/dx; Y_plot = 1 + 0.45/dy;
        T_plot(:,1) = T(X_plot,Y_plot,1);
        
        %Implementation of ADI
        T_inter = zeros(M-2,N-2); %Temperature at half intermediate time step
        T_store = zeros(M,N); %Store the obtained result of ADI
        for i = 1:nof_t-1
            
            Q1 = (dt/2)*Q(X,Y,(i-1/2)*dt); %Source term at next half time step
            Q2 = (dt/2)*Q(X,Y,i*dt);       %Source term at next time step
            s = 1;                  %Parameter used to execute two scheme turn by turn
            
            if s == 1
                b = B1*T(2:M-1,2:N-1,i) + Q1(2:M-1,2:N-1);
                %Since A1 act on the rows of T (x direction),  then the returned
                %soltion of this for loop is the transpose matrix of the actual one
                T_inter = thomas(A1,b);
                %B2 act on the rows of actual T_n+1/2, same effect as acting on
                %columns of T_inter from the last loop.
                b = B2*T_inter + Q2(2:M-1,2:N-1);
                T_inter = thomas(A2,b);
                %Return T_n+1
                T_store(2:M-1,2:N-1) = T_inter;
                s = s - 1;
                
            else %Reverse the scheme of x and y
                %F1 should act on the rows of T_n, the effect of which is same as
                %F1 acting on the columns of T_n'
                b = F1*T(2:M-1,2:N-1,i).' + Q1(2:M-1,2:N-1);
                T_inter = thomas(E1,b);
                %Return T_n+1/2
                b = F2*T_inter + Q2(floor(2:M-1.2:N-1));
                %E2 should act on the rows of T_n+1, the effect of which is same as
                %E2 acting on the columns of (T_n+1/2)'
                T_inter = thomas(E2,b);
                %Return (T_n+1)'
                T_inter = T_inter.';
                T_store(2:M-1,2:N-1) = T_inter;
                s = s + 1;
            end
            %Update Solution
            T(:,:,i+1) = T_store;
            %Solution at the next time step at x = 0.55, y = 0.45
            T_plot(10*(ac-1)+oc,i+1) = T(X_plot,Y_plot,i+1);
        end
        Plot(oc) = subplot(5,2,oc);  plot(t,T_plot(10*(ac-1)+oc,:)); grid on;
        title(['\alpha = ',num2str(alpha(ac)),' \Omega = ',num2str(O(oc))]);
        %Initialze the time vector and Temperature Matrix for the next loop
    end
    y_axes = get(gca,'YLim'); y_axes = 1.3*y_axes;
    ylim(Plot, y_axes);
    fig_index = sprintf('%d', ac);
    saveas(gcf,fig_index,'bmp')
end
toc

%Find out how alpha will affect the temperature revolution
P1 = figure; hold on ;
plot(t,T_plot(27,:),'k')
plot(t,T_plot(47,:),'--k')
plot(t,T_plot(97,:),'-.k')
xlabel('Time'); ylabel('Temperature');
title('Temperature Revolution at different \alpha with \Omega = 20.33');
legend('\alpha = 3.36','\alpha = 4.46','\alpha = 5');
y_axes = get(gca,'YLim'); ylim(1.1*y_axes);

%Next find out how omage will affect the temperature distribution
P2 = figure; hold on ;
plot(t,T_plot(12,:),'k')
plot(t,T_plot(14,:),'--k')
plot(t,T_plot(19,:),':k')
xlabel('Time'); ylabel('Temperature');
title('Temperature Revolution at different \Omega with \alpha = 0.67');
legend('\Omega = 4.22','\Omega = 10.67','\Omega = 30');
y_axes = get(gca,'YLim'); ylim(1.1*y_axes);



%Thomas Algorithmn
function solution = thomas(A,b)
s = size(A,1); %Number of rows
%Forward Elimination of subdiagonal elements;
for j = 1:s-1;
    
    M(j) = -A(j+1,j)/A(j,j);%Multiplier
    A(j+1,j) = A(j,j)*M(j) + A(j+1,j);
    A(j+1,j+1) = A(j+1,j+1)+A(j,j+1)*M(j);%New next diagonal term
    b(j+1,:) = b(j+1,:)+M(j)*b(j,:);
    %No need to modify subdiagonal terms since they are not used in the
    %following steps
    j=j+1;
end
%Backward substitution for solution;
solution(s,:) = b(s,:)/A(s,s); %The Last element of solution
for k = s-1:-1:1
    solution(k,:) = (b(k,:)-A(k,k+1)*solution(k+1,:))/A(k,k);%Update b vector with solutions  gradually
end
end