clear all, close all, clc
addpath('./utils/');
path2figs = '../Figures/ASYM_POT_WELL/'; mkdir(path2figs)
path2data = '../Data/'; mkdir(path2data)
ModelName = 'AsymPotentialWell_';

% Parameters
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
a = -0.25; 
B = [0; 0];
tspan = .001:.001:50;
x0 = [-1.50984, -1; 
      -1.1,      0; 
       a,    0.001];
Q  = 1;
R  = 1;
xREF = [1;0];

% Define functions
f = @(t,x,u)([x(2); -(x(1).^3 - x(1) - a.*x(1).^2 + a)]+B*u);
H = @(x)(x(2).^2/2 + x(1).^4/4 - x(1).^2/2 - (a/3).*x(1).^3 + a*x(1));
REF = H(xREF);

% Hamiltonian function
x = [-3:0.01:3]; 
y = [-3:0.01:3];
[X,Y] = meshgrid(x,y);
Hgrid = zeros(size(X));
for i = 1:size(X,1)
    for j = 1:size(X,2)
        Hgrid(i,j) = H([X(i,j); Y(i,j)]);
    end
end

%% Trajectory integration of unforced system for different initial
% conditions
y0 = zeros(length(tspan),2,size(x0,1));
Hvals0 = zeros(length(tspan),size(x0,1));
Jvals0 = zeros(length(tspan),size(x0,1));
for i = 1:size(x0,1)
    [t,y0tmp] = ode45(@(t,x)f(t,x,0),tspan,x0(i,:),ode_options);
    [Hvals0(:,i),Jvals0(:,i)] = evalCostFun_Hamiltonian(H,y0tmp,zeros(1,size(y0tmp,1)),Q,R,REF);
    y0(:,:,i) = y0tmp;
end

%% Controlled system // Parameters
gradH = @(x)([(x(1).^3 - x(1) - a.*x(1).^2 + a); x(2)]);

% B = [0; 1]; R = 1;
% ModelName1 = [ModelName, 'B01_'];
% B = [1; 0]; R= 1;
% ModelName1 = [ModelName, 'B10_'];
B = eye(2); R = 0.05.*eye(2);%0.05.*eye(2);
ModelName1 = [ModelName, 'B11_'];
Bc = B;

f       = @(t,x,u)([x(2); -(x(1).^3 - x(1) - a.*x(1).^2 + a)]+B*u);
gain    = @(x)(lqr(0,(gradH(x)'*B),Q,R));

%% 3-phase controller
[t,ytmp] = ode45(@(t,x)AsymmetricDoubleWellFUN(t,x,gain,H,a,B),tspan,x0(1,:),ode_options);
y1(:,:,1) = ytmp; 
[t,ytmp] = ode45(@(t,x)AsymmetricDoubleWellFUN(t,x,gain,H,a,B),tspan,x0(2,:),ode_options);
y1(:,:,2) = ytmp; 

%%
Hvals1 = zeros(length(t),size(y1,3));
Jvals1 = zeros(length(t),size(y1,3));
uvals1 = zeros(size(B,2),length(t),size(y1,3));
for i = 1:size(y1,3)
     [Hvals_tmp,Jvals_tmp,uvals_tmp] = evalCostFun_AsymPotentialWell_PhaseController(y1(:,:,i),B,R,Q,H,gain,a,REF);
     Hvals1(:,i) = Hvals_tmp;
     Jvals1(:,i) = Jvals_tmp;
     uvals1(:,:,i) = uvals_tmp;
end

%% SAVE RESULTS
save([path2data,[ModelName1,'Data.mat']])


