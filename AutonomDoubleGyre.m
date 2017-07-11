clear all, close all, clc

addpath('./utils/');
path2figs = '../Figures/DOUBLE_GYRE/'; mkdir(path2figs)
path2data = '../Data/'; mkdir(path2data)
ModelName = 'AutonomDoubleGyre_';

% Parameters double gyre
tspan   = 0:.01:50;
epsilon = 0;
omega   = 2*pi; 
A       = 0.25;

% Streamfunction / Hamiltonian
x = [0:0.01:2]; y = [0:0.01:1];
[X,Y]       = meshgrid(x,y);
StreamFun   = A*sin(pi*X).*sin(pi*Y);
Hfield      = -StreamFun;

% Parameters control
x0      = [1.3; 0.5];   % initial state
Q       = 1;            % weighing efun/state in cost function
R       = 1;            % penalize control
REF     = 0.2;         % reference value / level set
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);

% Define functions
Psi = @(x)(A*sin(pi*x(1)).*sin(pi*x(2)));           % Streamfunction
gradPsi = @(x)([A*pi*cos(pi*x(1))*sin(pi*x(2));     % Gradient of Psi
                A*pi*sin(pi*x(1))*cos(pi*x(2))]);

%% UNFORCED SYSTEM
B = [0; 0];
f = @(t,x,u)([-A*pi*sin(pi*x(1))*cos(pi*x(2));A*pi*cos(pi*x(1))*sin(pi*x(2))] + B*u);          
[t,y0] = ode45(@(t,x)f(t,x,0),tspan,x0,ode_options);
[Psivals0,Jvals0] = evalCostFun_KoopEfun(Psi,y0,zeros(1,size(y0,1)),Q,R,REF);

%% Controlled system: different actuation choice
% B = [0; 1]; R = 1;
% ModelName1 = [ModelName, 'B01_'];
% B = [1; 0]; R = 1;
% ModelName1 = [ModelName, 'B10_']; % not gain comp
B = [1, 0; 0, 1]; R  = eye(2);
ModelName1 = [ModelName, 'B11_'];
Bc = B;

% KRONIC using stream function
f = @(t,x,u)([-A*pi*sin(pi*x(1))*cos(pi*x(2)); A*pi*cos(pi*x(1))*sin(pi*x(2))] + B*u); % need to be updated with B
gain = @(x)(lqr(0,(gradPsi(x)'*B),Q,R));
[~,y1] = ode45(@(t,x)f(t,x,-gain(x)*(Psi(x)-REF)),tspan,x0,ode_options);
uvals1 = zeros(size(B,2),length(y1)); for k=1:length(y1), uvals1(:,k) = - gain(y1(k,:))*(Psi(y1(k,:))-REF); end
[Psivals1,Jvals1] = evalCostFun_KoopEfun(Psi,y1,uvals1,Q,R,REF);

%% Do the same for an ensemble of drifters
% using rk4 with constant timestep as it is faster

% Ensemble of initial conditions 
y0ic        = [0 0];
xvec        = 0.01:0.25:2; xvec = [xvec,1.98]; %0.5
yvec        = 0.01:0.25:1; yvec = [yvec,0.98]; %0.5
dt          = 0.005;
Ny          = length(yvec);
Nx          = length(xvec);
yIC         = zeros(2,Ny,Nx);
[x0,y0]     = meshgrid(xvec+y0ic(1),yvec+y0ic(2));
yIC(1,:,:)  = x0;
yIC(2,:,:)  = y0;

% UNFORCED SYSTEM
B       = [0; 0];
duration = 25;
tspan   =[0,duration]; 
L       = duration/dt;
yin     = yIC; 
yout    = zeros(2,Ny,Nx,L);
yout(:,:,:,1) = yin;
for step = 1:L-1
    time = step*dt; 
    yout(:,:,:,step+1) = rk4singlestep(@(t,y,u)DoubleGyre_ensemble(t,y,A,omega,epsilon,B,u),dt,time,yin,0);
    yin = yout(:,:,:,step);   
end

% CONTROL
% B = Bc;
% B = Bc;
B = Bc; Ncontrol = size(B,2);
gain = @(x)(lqr(0,(gradPsi(x)'*B),Q,R));

yin     = yIC; 
uin     = zeros(Ncontrol,Ny,Nx);
for iy = 1:Ny
    for ix = 1:Nx
        uin(:,iy,ix) = -gain(yin(:,iy,ix))*(Psi(yin(:,iy,ix))-REF);
    end
end
yout_ctrl = zeros(2,Ny,Nx,L);
yout_ctrl(:,:,:,1) = yin;
uout      = zeros(Ncontrol,Ny,Nx,L);
uout(:,:,:,1) = uin;
for step = 1:L-1
    time = step*dt; 
    yout_ctrl(:,:,:,step+1) = rk4singlestep(@(t,y,u)DoubleGyre_ensemble(t,y,A,omega,epsilon,B,u),dt,time,yin,uin);
    yin = yout_ctrl(:,:,:,step+1);  
    uin     = zeros(Ncontrol,Ny,Nx);
    for iy = 1:Ny
        for ix = 1:Nx
            uin(:,iy,ix) = -gain(yin(:,iy,ix))*(Psi(yin(:,iy,ix))-REF); 
        end
    end
    uout(:,:,:,step+1) = uin;
end

%% SAVE RESULTS
save([path2data,[ModelName1,'Ensemble.mat']])

