clear all , close all , clc 
path2figs = './../Figures/SLOW_MANIFOLD/'; mkdir(path2figs)

%% System
mu          = -0.1;
lambda      = -1;
A           = [mu 0 0; 0 lambda -lambda; 0 0 2*mu]; % Koopman linear dynamics
[T,D]       = eig(A);
slope_stab_man = T(3,3)/T(2,3); % slope of stable subspace (green)

%% Nonlinear phase plot with ensemble of trajectories // Case I
% STABLE

% Initial condition 
y0ic    = [0 0];
xvec    = -4:1:4;
yvec    = [-5,20];
x       = [-4:0.01:4];
dt      = 0.01;
Ny0     = length(yvec);
Nx0     = length(xvec);
yIC     = zeros(2,Ny0,Nx0);
[x0,y0] = meshgrid(xvec+y0ic(1),yvec+y0ic(2));
yIC(1,:,:) = x0;
yIC(2,:,:) = y0;

% System
mu          = -.1;
lambda      = -1;
duration    = 50; % for stable
tspan       = [0,duration]; 
L           = duration/dt;
yin         = yIC; 

yout = zeros(2,Ny0,Nx0,L);
yout(:,:,:,1) = yin;
for step = 1:L-1
    time = step*dt;
    yout(:,:,:,step+1) = rk4singlestep(@(t,y,u)SlowManifold_ensemble(t,y,mu,lambda),dt,time,yin,0);
    yin = yout(:,:,:,step);   
end

% UNSTABLE

% Initial condition
xvec    = -4:1:4;
Ny      = length(yvec);
Nx      = length(xvec);
yIC     = zeros(2,Ny,Nx);
[x0,y0] = meshgrid(xvec+y0ic(1),yvec+y0ic(2));
yIC(1,:,:) = x0;
yIC(2,:,:) = y0;

% Add control
mu          = -.1;
lambda      = 1;
ModelName   = 'SlowManifold_B01_';
B = [0; 1];
A = [mu 0; 0 lambda]; 
R = 1;
b = lambda/(lambda-2*mu);
T = [1 0 0; 0 1 -b; 0 0 1];
Aphi = [mu 0 0; 0 lambda 0; 0 0 2*mu];
Bphi = [B; 0];
Qphi = [1 0 0; 0 1 b; 0 b b^2];
Cphi = lqr(Aphi,Bphi,Qphi,R)*T;
vf3 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(Cphi(1:2)*x + Cphi(3)*x(1)^2);

SlowManifold_PhasePlot_Control

% Add control
mu          = -mu;
lambda      = -lambda;
ModelName   = [path2figs,'SlowManifold_B10_'];
subindex    = @(A,r,c) A(r,c);   
B = [1; 0];
A = [mu 0; 0 lambda]; 
R = 1;
Q = eye(2);
b = lambda/(lambda-2*mu);
T = [1 0 0; 0 1 -b];
Aphi = [mu 0; 0 lambda];
gradphi = @(x) ([1 0; -b*x(1) 1]);
gain = @(x)(lqr(Aphi,(gradphi(x)*B),Q,R)*T);
vf3 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(subindex(gain(x),1,1:2)*x + subindex(gain(x),1,3)*x(1)^2);

SlowManifold_PhasePlot_Control

%% Koopman eigenfunction optimal control
clear all , close all , clc
path2figs = './../Figures/SLOW_MANIFOLD/';
dt = 0.01;
tspan       = 0:dt:50; %50 
x0          = [-5; 5];
subindex    = @(A,r,c) A(r,c);      %# An anonymous function to index a matrix

%% OPTIMAL CONTROL B = [0; 1]
ModelName   = [path2figs, 'SlowManifold_B01_'];
B           = [0; 1];
mu          = -.1;
lambda      = 1;
b           = lambda/(lambda-2*mu);

% LQR on linearized system
A = [mu 0; 0 lambda]; 
Q = eye(2);
R = 1;
C = lqr(A,B,Q,R);
vf = @(t,x) A*x + [0; -lambda*x(1)^2] - B*C*x; 
[t,xLQR] = ode45(vf,tspan,x0);
uLQR = (C*xLQR')';
JLQR_x   = cumsum(xLQR(:,1).^2 + xLQR(:,2).^2 + uLQR.^2)'; 
JLQR_phi = cumsum(xLQR(:,1).^2 + xLQR(:,2).^2 + uLQR.^2 + ... % JLQR_x part
                 (xLQR(:,2)-b*xLQR(:,1).^2).^2); % via phi

% Koopman operator optimal control (KOOC); i.e., LQR on Koopman operator
A2 = [mu 0 0; 0 lambda -lambda; 0 0 2*mu]; 
B2 = [B; 0];
Q2 = [1 0 0; 0 1 0; 0 0 0];
R = 1;
C2 = lqr(A2,B2,Q2,R); % CONSTANT
% note that controller is nonlinear in the state x1
vf2 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(C2(1:2)*x + C2(3)*x(1)^2);
[t,xKOOC] = ode45(vf2,tspan,x0);
uKOOC = (C2(1:2)*xKOOC')'+C2(3)*xKOOC(:,1).^2;
JKOOC_x = cumsum(xKOOC(:,1).^2 + xKOOC(:,2).^2 + uKOOC.^2)'; 
JKOOC_phi = cumsum( xKOOC(:,1).^2 + xKOOC(:,2).^2 + uKOOC.^2 + ... % J_x part
                   (xKOOC(:,2)-b*xKOOC(:,1).^2).^2); % via phi

% Check
[Vy,Dy] = eig(A2');
Vy^(-1)*A2'*Vy

% "Analytic"
H2 = [A2 -B2*R^(-1)*B2';
     -Q2 -A2'];
[V2,D2] = eig(H2);
[D2,IX2] = sort(diag(D2),'ascend');
V2 = V2(:,IX2); D2 = D2(IX2);
P2 = V2(4:6,1:3)*V2(1:3,1:3)^(-1);
K2 = R^(-1)*B2'*P2

% KRONIC with eigenfunctions
b = lambda/(lambda-2*mu);
T = [1 0 0; 0 1 -b; 0 0 1]; % Transformation between y and phi
Aphi = [mu 0 0; 0 lambda 0; 0 0 2*mu];
Bphi = [B; 0];
Qphi = [1 0 0; 0 1 b; 0 b b^2];
%Qphi = [1 0 0; 0 1 0; 0 0 1];
Cphi = lqr(Aphi,Bphi,Qphi,R)*T; 
vf3 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(Cphi(1:2)*x + Cphi(3)*x(1)^2);
[t,xKOOC2] = ode45(vf3,tspan,x0);
uKOOC2 = (Cphi(1:2)*xKOOC2')' + Cphi(3)*xKOOC2(:,1).^2;
JKOOC2_x = cumsum(xKOOC2(:,1).^2 + xKOOC2(:,2).^2 + uKOOC2.^2 )';
JKOOC2_phi = cumsum(xKOOC2(:,1).^2 + xKOOC2(:,2).^2 + uKOOC2.^2 + ... % J_x part
                   (xKOOC2(:,2)-b*xKOOC2(:,1).^2).^2); % via phi
           
             
% "Analytic"
N = size(Aphi,1);
H3 = [Aphi -Bphi*R^(-1)*Bphi';
     -Qphi -Aphi'];
[V3,D3] = eig(H3);
[D3,IX3] = sort(diag(D3),'ascend');
V3 = V3(:,IX3); D3 = D3(IX3);
P3 = V3(N+1:2*N,1:N)*V3(1:N,1:N)^(-1);
K3 = R^(-1)*Bphi'*P3*T % equal to Cphi
% P*Aphi + Aphi'*P +Q - P*Bphi*R^(-1)*Bphi'*P

% Nonlinear Control
vf4 = @(t,x) [A*x(1:2) + [0; -lambda*x(1)^2] - B*R^(-1)*B'*x(3:4); ...
              -Q*x(1:2) - A'*x(3:4) - [0; lambda*x(1)]];
vf4_bc = @(xa,xb) [xa(1:2)-x0; xb(3:4)-eye(2)*xb(1:2)];
vf4_init = @(t) [x0; 0; 0]; 
solinit = bvpinit(tspan, vf4_init);
sol     = bvp4c(vf4,vf4_bc,solinit);
xNLC = sol.y(1:2,:)';
tNLC = sol.x;
uNLC = R^(-1)*B'*sol.y(3:4,:); uNLC = uNLC';
JNLC_x      = cumsum(xNLC(:,1).^2 + xNLC(:,2).^2 + (uNLC).^2)'; 
JNLC_phi    = cumsum(xNLC(:,1).^2 + xNLC(:,2).^2 + (uNLC).^2 + ... % J_x part
                     (xNLC(:,2)-b*xNLC(:,1).^2).^2);               % via phi
 
%A*sol.y(1:2,1) + [0; -lambda*sol.y(1,1)^2] - B*R^(-1)*B'*sol.y(3:4,1)


% Feedback linearization
A5 = [mu 0; 0 lambda]; 
Q5 = eye(2);
R5 = 1;
C5 = lqr(A5,B,Q5,R5);
vf5 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(-lambda*x(1)^2 + C5*x); 
[t,xFL] = ode45(vf5,tspan,x0);
uFL = (-lambda*(xFL(:,1).^2)' + C5*xFL')';
JFL_x   = cumsum(xFL(:,1).^2 + xFL(:,2).^2 + uFL.^2)'; 
JFL_phi = cumsum(xFL(:,1).^2 + xFL(:,2).^2 + uFL.^2 + ... % JLQR_x part
                (xFL(:,2)-b*xFL(:,1).^2).^2);               % via phi

% Plot (FIGURE 5)
showThree = 1;
ylim_vals = [-15 10];
axis_lim = [0 50 0 600000];
% axis_lim = [0 50 0 100000000];
showResults_SlowManifold

%% Balanced model
sys_y = ss(A2,B2,eye(3),0);
tf_y = tf(sys_y);
[sys_y_bal,g] = balreal(sys_y);
rsys = balred(sys_y,2);

sys_phi = ss(Aphi,Bphi,eye(3),0);
[sys_phi_bal,g_phi] = balreal(sys_phi);
rsys_phi = balred(sys_phi,2);

%% Cost function Q = eye (FIGURE 6)
% Koopman operator optimal control (KOOC)
A2 = [mu 0 0; 0 lambda -lambda; 0 0 2*mu]; 
B2 = [B; 0];
Q2 = [1 0 0; 0 1 0; 0 0 1];
R = 1;
C2 = lqr(A2,B2,Q2,R); % CONSTANT
vf2 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(C2(1:2)*x + C2(3)*x(1)^2);
[t,xKOOC] = ode45(vf2,tspan,x0);
uKOOC = (C2(1:2)*xKOOC')'+C2(3)*xKOOC(:,1).^2;
JKOOC_x = cumsum(xKOOC(:,1).^2 + xKOOC(:,2).^2 + uKOOC.^2)'; 
JKOOC_phi = cumsum( xKOOC(:,1).^2 + xKOOC(:,2).^2 + uKOOC.^2 + ... % J_x part
                 -2*b*xKOOC(:,1).^2.*xKOOC(:,2) + b^2*xKOOC(:,1).^4)'; % via phi
                 
% KC with eigenfunctions
b = lambda/(lambda-2*mu);
T = [1 0 0; 0 1 -b; 0 0 1]; % Transformation between y and phi
Aphi = [mu 0 0; 0 lambda 0; 0 0 2*mu];
Bphi = [B; 0];
Qphi = [1 0 0; 0 1 0; 0 0 1];
% T = [1 0 0; 0 1 -b];
% Aphi = [mu 0 ; 0 lambda];
% Bphi = [B];
% Qphi = [1 0; 0 1];
Cphi = lqr(Aphi,Bphi,Qphi,R)*T;
vf3 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(Cphi(1:2)*x + Cphi(3)*x(1)^2);
[t,xKOOC2] = ode45(vf3,tspan,x0);
uKOOC2 = (Cphi(1:2)*xKOOC2')' + Cphi(3)*xKOOC2(:,1).^2;
JKOOC2_x = cumsum(xKOOC2(:,1).^2 + xKOOC2(:,2).^2 + uKOOC2.^2 )';
JKOOC2_phi = cumsum(xKOOC2(:,1).^2 + xKOOC2(:,2).^2 + uKOOC2.^2 + ... % J_x part
                 - 2*b*xKOOC2(:,1).^2.*xKOOC2(:,2) + b^2*xKOOC2(:,1).^4)'; % via phi


% Show comparison
figure,hold on, box on
plot(tspan,JLQR_x,'-','Color',colors(1,:),'LineWidth',2)
plot(tspan,JKOOC_x,'-','Color',colors(2,:),'LineWidth',2)
plot(tspan,JKOOC2_x,'--','Color',colors(3,:),'LineWidth',2)
% plot(tspan,JNLC_x,':','Color','r')
xlabel('t'), ylabel('Jx')
axis([0 50 0 5*10^5])
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [ModelName,'Cost_Jx','.eps']);

figure,hold on, box on
plot(tspan,JLQR_phi,'-','Color',colors(1,:),'LineWidth',2)
plot(tspan,JKOOC_phi,'-','Color',colors(2,:),'LineWidth',2)
plot(tspan,JKOOC2_phi,'--','Color',colors(3,:),'LineWidth',2)
% plot(tspan,JNLC_phi,':','Color','r')
xlabel('t'), ylabel('Jp')
axis([0 50 0 8*10^5])
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [ModelName,'Cost_Jphi','.eps']);

%% Control
clear all , close all , clc
dt = 0.01;
tspan = 0:dt:50; 
x0 = [-5; 5];
subindex = @(A,r,c) A(r,c);      %# An anonymous function to index a matrix

path2figs = './../Figures/SLOW_MANIFOLD/';
%% OPTIMAL CONTROL B = [1; 0]
axis_lim = [0 50 0 15000];

ModelName = [path2figs,'SlowManifold_B10_'];
B       = [1; 0];
mu      = 0.1;
lambda  = -1;
b       = lambda/(lambda-2*mu);

% LQR on linearized system
A = [mu 0; 0 lambda]; 
Q = eye(2);
R = 1;
try
    C = lqr(A,B,Q,R);
    vf = @(t,x) A*x + [0; -lambda*x(1)^2] - B*C*x; 
    [t,xLQR] = ode45(vf,tspan,x0);
catch
    warning('Linearized system not stabilizable.')
    xLQR = zeros(length(tspan),2);
end
uLQR = (C*xLQR')';
JLQR_x   = cumsum(xLQR(:,1).^2 + xLQR(:,2).^2 + uLQR.^2)'; 
JLQR_phi = cumsum(xLQR(:,1).^2 + xLQR(:,2).^2 + uLQR.^2 ... % JLQR_x part
            - 2*b*xLQR(:,1).^2.*xLQR(:,2) + b^2*xLQR(:,1).^4)';  
        
% Control with truncated Koopman system; i.e., LQR on Koopman "matrix" (fails)
A2 = [mu 0 0; 0 lambda -lambda; 0 0 2*mu]; 
Q2 = [1 0 0; 0 1 0; 0 0 0];
R = 1;
try
    C2 = @(x) lqr(A2,[B; 2*x(1)],Q2,R); 
    vf2 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(subindex(C2(x),1,1:2)*x + subindex(C2(x),1,3)*x(1)^2);
    [t,xKOOC] = ode45(vf2,tspan,x0); % fails at some point
catch
    warning('Truncated Koopman system not stabilizable.')
    xKOOC = zeros(length(tspan),2);
    C2 = zeros(1,3);
end
uKOOC = (C2(1:2)*xKOOC')'+C2(3)*xKOOC(:,1).^2;
JKOOC_x = cumsum(xKOOC(:,1).^2 + xKOOC(:,2).^2 + uKOOC.^2)'; 
JKOOC_phi = cumsum( xKOOC(:,1).^2 + xKOOC(:,2).^2 + uKOOC.^2 + ... % J_x part
                 -2*b*xKOOC(:,1).^2.*xKOOC(:,2) + b^2*xKOOC(:,1).^4)'; % via phi
             
% KRONIC 
b = lambda/(lambda-2*mu);
T = [1 0 0; 0 1 -b];
Aphi = [mu 0; 0 lambda];
gradphi = @(x) ([1 0; -2*b*x(1) 1]); 
gain = @(x)(lqr(Aphi,(gradphi(x)*B),Q,R)*T);
vf3 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(subindex(gain(x),1,1:2)*x + subindex(gain(x),1,3)*x(1)^2);
[t,xKOOC2] = ode45(vf3,tspan,x0);

uKOOC2 = zeros(length(t),1);
for i = 1:length(t)
    uKOOC2(i) = subindex(gain(xKOOC2(i,:)),1,1:2)*xKOOC2(i,:)' + subindex(gain(xKOOC2(i,:)),1,3)*xKOOC2(i,1)^2;
end
JKOOC2_x = cumsum(xKOOC2(:,1).^2 + xKOOC2(:,2).^2 + uKOOC2.^2 )';
JKOOC2_phi = cumsum(xKOOC2(:,1).^2 + xKOOC2(:,2).^2 + uKOOC2.^2 + ... % J_x part
                 - 2*b*xKOOC2(:,1).^2.*xKOOC2(:,2) + b^2*xKOOC2(:,1).^4)'; % via phi
             

% Nonlinear Control
vf4 = @(t,x) [A*x(1:2) + [0; -lambda*x(1)^2] - B*R^(-1)*B'*x(3:4); ...
              -Q*x(1:2) - A'*x(3:4) - [0; lambda*x(1)]];
vf4_bc = @(xa,xb) [xa(1:2)-x0; xb(3:4)-eye(2)*xb(1:2)];
vf4_init = @(t) [x0; 0; 0]; 
solinit = bvpinit(tspan, vf4_init);
sol     = bvp4c(vf4,vf4_bc,solinit);
xNLC = sol.y(1:2,:)';
tNLC = sol.x;
uNLC = R^(-1)*B'*sol.y(3:4,:); uNLC = uNLC';
JNLC_x = cumsum(xNLC(:,1).^2 + xNLC(:,2).^2 + (uNLC).^2)'; 
JNLC_phi = cumsum( xNLC(:,1).^2 + xNLC(:,2).^2 + (uNLC).^2 + ... % J_x part
                 -2*b*xNLC(:,1).^2.*xNLC(:,2) + b^2*xNLC(:,1).^4)'; % via phi

% Feedback linearization (fails, set to zero)
% A5 = [0 1; -2*mu*lambda lambda+2*mu];
% B5 = [0;1];
% Q5 = eye(2);
% R5 = 1;
% C5 = lqr(A5,B5,Q5,R5);
% vf5 = @(t,x) A*x + [0; -lambda*x(1)^2] + B*(1/(2*lambda*x(1))*C5*x); 
% [t,xFL] = ode45(vf5,tspan,x0);             
xFL = zeros(length(tspan),2);
C5 = zeros(1,2);
uFL = ((1./(2*lambda*xFL(:,1)')).*(C5*xFL'))';
JFL_x   = cumsum(xFL(:,1).^2 + xFL(:,2).^2 + uFL.^2)'; 
JFL_phi = cumsum(xFL(:,1).^2 + xFL(:,2).^2 + uFL.^2 ... % JLQR_x part
            - 2*b*xFL(:,1).^2.*xFL(:,2) + b^2*xFL(:,1).^4)';  

ylim_vals = [-5 10];
showThree = 1;
showResults_SlowManifold

