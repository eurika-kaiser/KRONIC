clear all, close all, clc

ModelName = 'Pendulum_';
path2data = '../Data/'; mkdir(path2data)
path2figs = '../Figures/PENDULUM/'; mkdir(path2figs)

% Parameters
tspan = .001:.001:100;
x0 = [0.5;0.5];[.5; .5];
Q  = 1;
R  = 1;
H = @(x)([(1/2)*x(2)^2-cos(x(1))]);
gradH = @(x)([sin(x(1)); x(2)]);
[X,Y] = meshgrid([-2*pi:0.01:2*pi], [-4:0.01:4]);
Hfield = (1/2)*Y.^2-cos(X);

%% Unforced
REF = 0;
B = [0; 0];
f = @(t,x,u)([x(2); -sin(x(1))]+B*u);
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
[t,y0] = ode45(@(t,x)f(t,x,0),tspan,x0,ode_options);
[Hvals0,Jvals0] = evalCostFun_Hamiltonian(H,y0,zeros(1,size(y0,1)),Q,R,REF);

%% Controlled system
% B = [0; 1];
% ModelName1 = [ModelName, 'B01_'];
B = [1; 0];
ModelName1 = [ModelName, 'B10_'];
f = @(t,x,u)([x(2); -sin(x(1))]+B*u);


% KRONIC: Case 1
REF = 1;%H([pi;0]);
Href = @(x)([(1/2)*x(2)^2-cos(x(1))]-REF);
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[~,y1] = ode45(@(t,x)f(t,x,-gain(x)*Href(x)),tspan,x0,ode_options);
uvals1 = zeros(1,length(y1)); for k=1:length(y1), uvals1(1,k) = - gain(y1(k,:))*Href(y1(k,:)); end
[Hvals1,Jvals1] = evalCostFun_Hamiltonian(H,y1,uvals1,Q,R,REF);

% Store results
DataStore.y1 = y1;
DataStore.u1 = uvals1;
DataStore.H1 = Hvals1;
DataStore.J1 = Jvals1;
DataStore.tspan1 = tspan;


% KRONIC: Case 2
REF = 0;
Href = @(x)([(1/2)*x(2)^2-cos(x(1))]-REF);
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[~,y1] = ode45(@(t,x)f(t,x,-gain(x)*Href(x)),tspan,x0,ode_options);
uvals1 = zeros(1,length(y1)); for k=1:length(y1), uvals1(1,k) = - gain(y1(k,:))*Href(y1(k,:)); end
[Hvals1,Jvals1] = evalCostFun_Hamiltonian(H,y1,uvals1,Q,R,REF);

% Store results
DataStore.y2 = y1;
DataStore.u2 = uvals1;
DataStore.H2 = Hvals1;
DataStore.J2 = Jvals1;
DataStore.tspan2 = tspan;

% KRONIC: Case 3
REF = 2;
Href = @(x)([(1/2)*x(2)^2-cos(x(1))]-REF);
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[~,y1] = ode45(@(t,x)f(t,x,-gain(x)*Href(x)),tspan,x0,ode_options);
uvals1 = zeros(1,length(y1)); for k=1:length(y1), uvals1(1,k) = - gain(y1(k,:))*Href(y1(k,:)); end
[Hvals1,Jvals1] = evalCostFun_Hamiltonian(H,y1,uvals1,Q,R,REF);

Htest = zeros(size(Hvals1));
for k=1:length(y1), Htest(k) = H([mod(y1(k,1),2*pi), y1(k,2)]); end % seems ok, =2

% Store results
DataStore.y3 = y1;
DataStore.u3 = uvals1;
DataStore.H3 = Hvals1;
DataStore.J3 = Jvals1;
DataStore.tspan3 = tspan;


% KRONIC: Case 4
REF = H([0,0]);
Href = @(x)([(1/2)*x(2)^2-cos(x(1))]-REF);
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[~,y1] = ode45(@(t,x)f(t,x,-gain(x)*Href(x)),tspan,x0,ode_options);
uvals1 = zeros(1,length(y1)); for k=1:length(y1), uvals1(1,k) = - gain(y1(k,:))*Href(y1(k,:)); end
[Hvals1,Jvals1] = evalCostFun_Hamiltonian(H,y1,uvals1,Q,R,REF);

% Store results
DataStore.y4 = y1;
DataStore.u4 = uvals1;
DataStore.H4 = Hvals1;
DataStore.J4 = Jvals1;
DataStore.tspan4 = tspan;

%% SAVE RESULTS
save([path2data,[ModelName1,'Data.mat']])

