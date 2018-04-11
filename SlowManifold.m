clear all , close all , clc 
path2figs = './../Figures/SLOW_MANIFOLD/'; mkdir(path2figs)
addpath('./utils');

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
ModelName   = [path2figs,ModelName];
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
tspan       = 0:dt:50; 
x0          = [-5; 5];
subindex    = @(A,r,c) A(r,c);      %# An anonymous function to index a matrix

%% OPTIMAL CONTROL B = [0; 1]
ModelName   = [path2figs, 'SlowManifold_B01_'];
B           = [0; 1];
mu          = -.1;
lambda      = 1;
b           = lambda/(lambda-2*mu);

% LQR on linearized system
A = [mu 0; 0 lambda];   % system matrix from linearization
Q = eye(2);             % state penalization
R = 1;                  % control penalization
C = lqr(A,B,Q,R);       % gain
vf = @(t,x) A*x + [0; -lambda*x(1)^2] - B*C*x; % apply gain in real system
[t,xLQR] = ode45(vf,tspan,x0); 
uLQR = (C*xLQR')';      % recalculate control input
JLQR_x   = cumsum(xLQR(:,1).^2 + xLQR(:,2).^2 + uLQR.^2)'; % standard LQR cost
% JLQR_phi = cumsum(xLQR(:,1).^2 + xLQR(:,2).^2 + uLQR.^2 + ... % cost in eigenfun, JLQR_x part
%                  (xLQR(:,2)-b*xLQR(:,1).^2).^2); % via phi
JLQR_phi = cumsum(xLQR(:,1).^2 + xLQR(:,1).^2 + uLQR.^2 + ... % cost in eigenfun,
                 (xLQR(:,2)-b*xLQR(:,1).^2).^2); % via phi             

             
% Koopman operator optimal control (KOOC); i.e., LQR on Koopman operator
A2 = [mu 0 0; 0 lambda -lambda; 0 0 2*mu]; % system matrix in (nonlinear) observables
B2 = [B; 0];                               % new control vector in observables
Q2 = [1 0 0; 0 1 0; 0 0 0];                % weight matrix only penalizing x1,x2
R = 1;
C2 = lqr(A2,B2,Q2,R); 
% note that controller becomes nonlinear in the state x1
vf2 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(C2(1:2)*x + C2(3)*x(1)^2);
[t,xKOOC] = ode45(vf2,tspan,x0);
uKOOC = (C2(1:2)*xKOOC')'+C2(3)*xKOOC(:,1).^2;
JKOOC_x = cumsum(xKOOC(:,1).^2 + xKOOC(:,2).^2 + uKOOC.^2)'; 
JKOOC_phi = cumsum( xKOOC(:,1).^2 + (xKOOC(:,1).^2).^2 + uKOOC.^2 + ... % J_x part
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
JKOOC2_phi = cumsum(xKOOC2(:,1).^2 + (xKOOC2(:,1).^2).^2 + uKOOC2.^2 + ... % J_x part
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
JNLC_phi    = cumsum(xNLC(:,1).^2 + (xNLC(:,1).^2).^2 + (uNLC).^2 + ... % J_x part
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
JFL_phi = cumsum(xFL(:,1).^2 + (xFL(:,1).^2).^2 + uFL.^2 + ... % JLQR_x part
                (xFL(:,2)-b*xFL(:,1).^2).^2);               % via phi

% Plot figures
showThree = 1;
ylim_vals = [-15 10];
axis_lim = [0 50 0 600000];
showResults_SlowManifold


%% Balanced model
% sys_y = ss(A2,B2,eye(3),0);
% tf_y = tf(sys_y);
% [sys_y_bal,g] = balreal(sys_y);
% rsys = balred(sys_y,2);
% 
% sys_phi = ss(Aphi,Bphi,eye(3),0);
% [sys_phi_bal,g_phi] = balreal(sys_phi);
% rsys_phi = balred(sys_phi,2);
% 
% [sysr,u] = minreal(sys_phi);

% LQR on minreal
% Q = eye(3);
% C = lqr(u*sys_phi.A*u',u*sys_phi.B,Q,R);       % gain
% vf = @(t,x) A*x + [0; -lambda*x(1)^2] - B*C*x; % apply gain in real system
% [t,xKD] = ode45(vf,tspan,x0); 

%% Cost function Q = eye 
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
JKOOC_phi = cumsum( xKOOC(:,1).^2 + (xKOOC(:,1).^2).^2 + uKOOC.^2 + ... % J_x part
                    (xKOOC(:,2)-b*xKOOC(:,1).^2).^2); % via phi
                
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
JKOOC2_phi = cumsum(xKOOC2(:,1).^2 + (xKOOC2(:,1).^2).^2 + uKOOC2.^2 + ... % J_x part
                    (xKOOC2(:,2)-b*xKOOC2(:,1).^2).^2); % via phi

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
tspan = 0:dt:50;%50; 
x0 = [-5; 5];
subindex = @(A,r,c) A(r,c);      %# An anonymous function to index a matrix

path2figs = './../Figures/SLOW_MANIFOLD/';
%% OPTIMAL CONTROL B = [1; 0]
axis_lim = [0 50 0 20000];

ModelName = [path2figs,'SlowManifold_B10_'];
B       = [1; 0];
mu      = 0.1;
lambda  = -1;
b       = lambda/(lambda-2*mu);
R       = 1;
RKRONIC = 4.0;
phifun3 = @(x) [x(1); x(2)-b*x(1)^2; x(1)^2];

% LQR on linearized system
A = [mu 0; 0 lambda];
Q = eye(2);
C = lqr(A,B,Q,R);
vf = @(t,x) A*x + [0; -lambda*x(1)^2] - B*C*x;
[t,xLQR] = ode45(vf,tspan,x0);
uLQR = (C*xLQR')';
phiLQR = zeros(length(t),3);
for i = 1:length(t)
    phiLQR(i,:) = phifun3(xLQR(i,:));
end
JLQR_x   = cumsum(xLQR(:,1).^2 + xLQR(:,2).^2 + R*uLQR.^2)';
JLQR_phi = cumsum(phiLQR(:,1).^2 + phiLQR(:,2).^2 + RKRONIC*uLQR.^2);
              
% Control with truncated Koopman system; i.e., LQR on Koopman "matrix" (fails)
A2 = [mu 0 0; 0 lambda -lambda; 0 0 2*mu];
Q2 = [1 0 0; 0 1 0; 0 0 0];

% Minreal
% B2 = [B; 2*x0(1)];
% sys_y = ss(A2,B2,eye(3),0);
% [sys_balreal,g,T,Ti] = balreal(sys_y);
% sys_balred = balred(sys_y,2);
% [sys_minreal,u_mr] = minreal(sys_y);
% tf_y = tf(sys_y);

try
    C2 = @(x) lqr(A2,[B; 2*x(1)],Q2,0.01*R);
    vf2 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(subindex(C2(x),1,1:2)*x + subindex(C2(x),1,3)*x(1)^2);
    [t,xKOOC] = ode45(vf2,tspan,x0); % fails at some point
    
%         B2 = @(x) [[B; 2*x(1)], [0;1;0]];
%         C2 = @(x) lqr(A2,B2(x),Q2,0.05*eye(2));
%         vf2 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(subindex(C2(x),1,1:2)*x + subindex(C2(x),1,3)*x(1)^2);
%         [t,xKOOC] = ode45(vf2,tspan,x0); % fails at some point
%         
%         B2 = @(x) [[B; 2*x(1)], [0;1;0]];
%         C2 = @(x) lqr(A2,B2(x),Q2,eye(2)); % getBRSys(sys_y,B2(x)),R,
%         T  = @(x) getBRTrafo(sys_y,B2(x));
%         vf2 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*C2(x)*T(x)*[x;x(1).^2];
%         [t,xKOOC] = ode45(vf2,tspan,x0); % fails at some point
    uKOOC = zeros(length(t),1);
    for i = 1:length(t)
        uKOOC(i) = subindex(C2(xKOOC(i,:)),1,1:2)*xKOOC(i,:)' + subindex(C2(xKOOC(i,:)),1,3)*xKOOC(i,1)^2;
    end
catch
    warning('Truncated Koopman system not stabilizable.')
    xKOOC = zeros(length(tspan),2);
    C2 = zeros(1,3);
    uKOOC = (C2(1:2)*xKOOC')'+C2(3)*xKOOC(:,1).^2; 
end

phiKOOC = zeros(length(t),3);
for i = 1:length(t)
    phiKOOC(i,:) = phifun3(xKOOC(i,:));
end
JKOOC_x = cumsum(xKOOC(:,1).^2 + xKOOC(:,2).^2 + R*uKOOC.^2)';
JKOOC_phi = cumsum(phiKOOC(:,1).^2 + phiKOOC(:,2).^2 + RKRONIC*uKOOC.^2);

% KRONIC 
% Truncated 2D
b = lambda/(lambda-2*mu);
T = [1 0 0; 0 1 -b]; % y -> phi
Aphi = [mu 0; 0 lambda];
gradphi = @(x) ([1 0; -2*b*x(1) 1]);
Qphi = eye(2);
gain = @(x)(lqr(Aphi,(gradphi(x)*B),Qphi,RKRONIC)*T);
vf3 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(subindex(gain(x),1,1:2)*x + subindex(gain(x),1,3)*x(1)^2);
% gain = @(x)(lqr(Aphi,(gradphi(x)*B),Qphi,R));
% vf3 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*gain(x)*phifun2(x); % above are identical to this
[t,xKOOC2] = ode45(vf3,tspan,x0);

% 3D
% b = lambda/(lambda-2*mu);
% T = [1 0 0; 0 1 -b; 0, 0, 1]; 
% Aphi = [mu 0 0; 0 lambda 0; 0 0 2*mu];
% gradphi = @(x) ([1 0; -2*b*x(1) 1; 2*x(1) 0]);
% % Qphi = eye(3);
% Qphi = [1 0 0; 0 1 b; 0 b b^2];
% gain = @(x)(lqr(Aphi,(gradphi(x)*B),Qphi,0.01*RKRONIC)*T); %RKRONIC
% vf3 = @(t,x) A*x + [0; -lambda*x(1)^2] - B*(subindex(gain(x),1,1:2)*x + subindex(gain(x),1,3)*x(1)^2);
% [t,xKOOC2] = ode45(vf3,tspan(1:10),x0);


uKOOC2 = zeros(length(t),1);
phiKOOC2 = zeros(length(t),3);
for i = 1:length(t)
    uKOOC2(i) = subindex(gain(xKOOC2(i,:)),1,1:2)*xKOOC2(i,:)' + subindex(gain(xKOOC2(i,:)),1,3)*xKOOC2(i,1)^2;
%     uKOOC2(i) = gain(xKOOC2(i,:))*phifun2(xKOOC2(i,:));
    phiKOOC2(i,:) = phifun3(xKOOC2(i,:));
end
JKOOC2_x = cumsum(xKOOC2(:,1).^2 + xKOOC2(:,2).^2 + R*uKOOC2.^2 )';
JKOOC2_phi = cumsum(phiKOOC2(:,1).^2 + phiKOOC2(:,2).^2 + RKRONIC*uKOOC2.^2);

                
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
phiNLC = zeros(length(t),3);
for i = 1:length(t)
    phiNLC(i,:) = phifun3(xNLC(i,:));
end
JNLC_x = cumsum(xNLC(:,1).^2 + xNLC(:,2).^2 + R*(uNLC).^2)';
JNLC_phi = cumsum(phiNLC(:,1).^2 + phiNLC(:,2).^2 + RKRONIC*uNLC.^2);
              
% Feedback linearization (fails, set to zero)
% A5 = [0 1; -2*mu*lambda lambda+2*mu];
% B5 = [0;1];
% Q5 = eye(2);
% R5 = 1;
% C5 = lqr(A5,B5,Q5,R5);
% vf5 = @(t,x) A*x + [0; -lambda*x(1)^2] + B*(1/(2*lambda*x(1))*C5*x); 
% [t,xFL] = ode45(vf5,tspan,x0);             
xFL = zeros(length(tspan),2); % feedback linear. doesn't work, set to zero for plotting
C5 = zeros(1,2);
uFL = ((1./(2*lambda*xFL(:,1)')).*(C5*xFL'))';
JFL_x   = cumsum(xFL(:,1).^2 + xFL(:,2).^2 + R*uFL.^2)'; 
phiFL = zeros(length(t),3);
for i = 1:length(t)
    phiFL(i,:) = phifun3(phiFL(i,:));
end
JFL_phi = cumsum(phiFL(:,1).^2 + phiFL(:,2).^2 + RKRONIC*uFL.^2);
            
ylim_vals = [-5 10];
showThree = 1;
showResults_SlowManifold

return
%% Different plots
figure, box on
plot(xNLC(:,1),xNLC(:,2),'-k','LineWidth',2), hold on
plot(xLQR(:,1),xLQR(:,2),'--r','LineWidth',2)
% plot(xKOOC(:,1),xKOOC(:,2),'-.m','LineWidth',2)
plot(xKOOC2(:,1),xKOOC2(:,2),':b','LineWidth',2)
title('Phase Plot')

figure, box on
plot(xNLC(:,1),'-k','LineWidth',2), hold on
plot(xNLC(:,2),'-k','LineWidth',2), 
plot(xLQR(:,1),'--r','LineWidth',2)
plot(xLQR(:,2),'--r','LineWidth',2)
plot(xKOOC2(:,1),':b','LineWidth',2)
plot(xKOOC2(:,2),':b','LineWidth',2)
% plot(xKOOC(:,1),'-.m','LineWidth',2)
% plot(xKOOC(:,2),'-.m','LineWidth',2)
title('Time series')

figure, box on
semilogy((JNLC_x),'-k','LineWidth',2), hold on
semilogy((JLQR_x),'--r','LineWidth',2)
semilogy((JKOOC2_x),':b','LineWidth',2)
semilogy((JKOOC_x),'-.m','LineWidth',2)
title('Jx')

figure, box on
semilogy((JNLC_phi),'-k','LineWidth',2), hold on
semilogy((JLQR_phi),'--r','LineWidth',2)
semilogy((JKOOC2_phi),':b','LineWidth',2)
semilogy((JKOOC_phi),'-.m','LineWidth',2)
title('Jphi')

figure, hold on, box on
plot(cumsum((uNLC).^2),'-k','LineWidth',2)
plot(cumsum((uLQR).^2),'--r','LineWidth',2)
plot(cumsum((uKOOC2).^2),':b','LineWidth',2)
plot(cumsum((uKOOC).^2),'-.m','LineWidth',2)
title('sum u^2')

% figure, hold on, box on
% plot(cumsum((uNLC).^2),'-k','LineWidth',2)
% plot(cumsum((uLQR).^2),'--r','LineWidth',2)
% plot(cumsum((uKOOC2).^2),':b','LineWidth',2)
%%
figure,
subplot(3,1,1)
hold on, box on, title('LQR')
plot(cumsum((uLQR).^2),'-k','LineWidth',2)
plot(cumsum(xLQR(:,1).^2),'--r','LineWidth',2)
plot(cumsum((xLQR(:,2)-b*xLQR(:,1).^2).^2),':b','LineWidth',2)
plot(JLQR_phi,'-g','LineWidth',2)
ylim([0 15000])
subplot(3,1,2)
hold on, box on, title('KRONIC')
plot(cumsum(RKRONIC*(uKOOC2).^2),'-k','LineWidth',2)
plot(cumsum(xKOOC2(:,1).^2),'--r','LineWidth',2)
plot(cumsum((xKOOC2(:,2)-b*xKOOC2(:,1).^2).^2),':b','LineWidth',2)
plot(JKOOC2_phi,'-g','LineWidth',2)
ylim([0 15000])
subplot(3,1,3)
hold on, box on, title('NLC')
plot(cumsum((uNLC).^2),'-k','LineWidth',2)
plot(cumsum(xNLC(:,1).^2),'--r','LineWidth',2)
plot(cumsum((xNLC(:,2)-b*xNLC(:,1).^2).^2),':b','LineWidth',2)
plot(JNLC_phi,'-g','LineWidth',2)
legend('u','x1','x2-bx1^2','J')
ylim([0 15000])

%%
xspace = [-5:0.1:5]; Nx = length(xspace);
[X,Y] = meshgrid(xspace,xspace);
PHI = zeros(Nx,Nx,3);
for ix = 1:Nx
    for iy = 1:Nx
        PHI(iy,ix,1:3) = phifun3([X(iy,ix); Y(iy,ix)]);
    end
end
figure,
subplot(1,3,1)
imagesc(xspace,xspace,PHI(:,:,1))
title('phi1')
subplot(1,3,2)
imagesc(xspace,xspace,PHI(:,:,2))
title('phi2')
subplot(1,3,3)
imagesc(xspace,xspace,PHI(:,:,3))
title('phi3')


