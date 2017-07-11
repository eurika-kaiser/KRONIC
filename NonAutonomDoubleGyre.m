clear all, close all, clc

path2figs = '../Figures/DOUBLE_GYRE/'; mkdir(path2figs)
path2data = '../Data/'; mkdir(path2data)
ModelName = 'NonAutonomDoubleGyre_';

% Parameters double gyre
tspan   = 0:.01:50;
epsilon = 0.25;
omega   = 2*pi;
A       = 0.25;

% Forcing in x
% x = 1+ (sqrt(1+4*epsilon^2*sin(omega*tspan).^2)-1)./(2*epsilon*sin(omega*tspan));
% figure,plot(tspan,x)

% Streamfunction / Hamiltonian
PSIsteady     = @(x)(A*sin(pi*x(1)).*sin(pi*x(2)));
Hsteady     = @(x)(-A*sin(pi*x(1)).*sin(pi*x(2)));
x = [0:0.01:2]; y = [0:0.01:1];
[X,Y]   = meshgrid(x,y);
force   = zeros(size(X,1),size(X,2),length(tspan));
PSIfield  = zeros(size(X,1),size(X,2),length(tspan));
for it = 1:length(tspan)
    force(:,:,it)       = epsilon*sin(omega.*tspan(it))*X.^2 + (1-2*epsilon*sin(omega.*tspan(it))).*X;
    PSIfield(:,:,it)    = A*sin(pi*force(:,:,it)).*sin(pi*Y);
end
StreamFun   = A*sin(pi*X).*sin(pi*Y);

% Parameters control
% x0 = [1.3; 0.5];
x0 = [1.02; 0.01];
Q  = 1;
R  = 1;
REF = 0.2;
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
% ode_options = odeset('RelTol',1e-6, 'AbsTol',1e-7);

% Define functions
force   = @(x,t)(epsilon*sin(omega*t)*x(1).^2 + (1-2*epsilon*sin(omega*t)).*x(1));
dfdx    = @(x,t)(2*epsilon*sin(omega*t)*x(1) + (1-2*epsilon*sin(omega*t)));
Psi     = @(x,t)(A*sin(pi*force(x,t)).*sin(pi*x(2)));
gradPsi = @(x,t)([A*pi*cos(pi*force(x,t))*sin(pi*x(2))*dfdx(x,t); ...
                  A*pi*sin(pi*force(x,t))*cos(pi*x(2))]);


%% UNFORCED SYSTEM
B = [0; 0];
f = @(t,x,u)([-A*pi*sin(pi*force(x,t))*cos(pi*x(2)); ...
    A*pi*cos(pi*force(x,t))*sin(pi*x(2)).*dfdx(x,t)]+B*u);  %TODO, check again, dfdx
[t,y0] = ode45(@(t,x)f(t,x,0),tspan,x0,ode_options);
[Psivals0,Jvals0] = evalCostFun_TimeDependentKoopEfun(Psi,y0,zeros(size(y0,1),1),Q,R,REF,tspan);

%% KRONIC using stream function
% B = [0; 1]; R = 1;
% ModelName = [ModelName, 'B01_'];
% B = [1; 0]; R = 1;
% ModelName = [ModelName, 'B10_'];
B = eye(2); R = eye(2);
ModelName1 = [ModelName, 'B11_'];
Bc = B;

% Koopman eigenfunction control
f = @(t,x,u)([-A*pi*sin(pi*force(x,t))*cos(pi*x(2)); ...
               A*pi*cos(pi*force(x,t))*sin(pi*x(2))]+B*u);
gain = @(x,t)(lqr(0,(gradPsi(x,t)'*B),Q,R));
[~,y1] = ode45(@(t,x)f(t,x,-gain(x,t)*(Psi(x,t)-REF)),tspan,x0,ode_options);

uvals1 = zeros(length(y1),size(B,2)); for k=1:length(y1), uvals1(k,:) = - gain(y1(k,:),tspan(k))*(Psi(y1(k,:),tspan(k))-REF); end
[Psivals1,Jvals1] = evalCostFun_TimeDependentKoopEfun(Psi,y1,uvals1,Q,R,REF,tspan);


%% SAVE RESULTS
save([path2data,[ModelName1,'Data.mat']])