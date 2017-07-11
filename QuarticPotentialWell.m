clear all, close all, clc
figpath = '../Figures/QUARTIC_POT_WELL/'; mkdir(figpath)
path2data = '../Data/'; mkdir(path2data)
addpath('./utils');

% Parameters
n = 2;
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,n));
ModelName = 'QuarticPotentialWell';
dt = .025;
tspan = 0:dt:10;

% Define functions
f = @(t,x,u)([x(2); -x(1)^3]);
H = @(x)((1/2)*x(2)^2 + (1/4)*x(1)^4);

%% Hamiltonian Fuction Schematic
[t,y1] = ode45(f,tspan,[-2 0],options); H1 = H([-2 0]);
[t,y2] = ode45(f,tspan,[-3 0],options); H2 = H([-3 0]);
[t,y3] = ode45(f,tspan,[-4 0],options); H3 = H([-4 0]);

y0 = [-4; -10];
[t,yb] = ode45(f,tspan,y0,options);
Hb = H(y0);

x = -5:.02:5; Nx = length(x);
y = -20:.02:20; Ny = length(y);
[X,Y] = meshgrid(x,y);
Hmat = (1/2)*Y.^2 + (1/4)*X.^4;
Hmat(Hmat(:)>Hb) = NaN;
Hmat = reshape(Hmat,[Nx*Ny, 1]);
XY = [reshape(X, [Nx*Ny, 1]), reshape(Y, [Nx*Ny, 1])];
TF = isnan(Hmat);
XY(TF,:) = [];

tri = delaunay(XY(:,1),XY(:,2));
Hxy = (1/2)*XY(:,2).^2 + (1/4)*XY(:,1).^4;

figure,
trihandle = trisurf(tri, XY(:,1), XY(:,2),Hxy); hold on, shading interp
colormap([0 0 0.7])
% surf(x,y,Hmat), shading interp, hold on
plot3(yb(:,1),yb(:,2),Hb.*ones(size(yb)),'-k','LineWidth',1)
plot3(y1(:,1),y1(:,2),H1.*ones(size(yb)),'-k','LineWidth',1)
plot3(y2(:,1),y2(:,2),H2.*ones(size(yb)),'-k','LineWidth',1)
plot3(y3(:,1),y3(:,2),H3.*ones(size(yb)),'-k','LineWidth',1)
plot3(0,0,H([0,0]),'ok','MarkerSize',5,'MarkerFaceColor','k')
% plot3(y1(:,1),y1(:,2),0.*ones(size(yb)),'-k','LineWidth',1)
% plot3(y2(:,1),y2(:,2),0.*ones(size(yb)),'-k','LineWidth',1)
% plot3(y3(:,1),y3(:,2),0.*ones(size(yb)),'-k','LineWidth',1)
xlabel('x')
ylabel('y')
zlabel('z')
trihandle.FaceAlpha = 0.4;
axis off
set(gcf,'Position',[100 100 400 286]);
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [figpath, ModelName,'_Energy.eps']);


%% Unforced system
ModelName = 'QuarticPotentialWell_REF002_';

dt = 0.01;
tspan = .001:dt:100;
Q  = 1;
R  = eye(2);
REF = 0.15; 
Ni = size(R,2);

% Unforced
B = zeros(2);
f = @(t,x,u)([x(2); -x(1)^3]+B*u);
H = @(x)((1/2)*x(2)^2 + (1/4)*x(1)^4 - REF);
% y = +/- sqrt(2*REF- 1/2 *x^4)
gradH = @(x)([x(1)^3; x(2)]);

%% KRONIC
B = eye(2);
ModelName1 = [ModelName, 'B11_'];
f = @(t,x,u)([x(2); -x(1)^3] + B*u);

% Koopman eigenfunction control
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[~,y1] = ode45(@(t,x)f(t,x,-gain(x)*H(x)),tspan,[.001; .001],ode_options);
gain = @(x)(lqr(0,(gradH(x)'*B),Q,R));
[~,y2] = ode45(@(t,x)f(t,x,-gain(x)*H(x)),tspan,[1; 1],ode_options);

%% SAVE RESULTS
save([path2data,[ModelName1,'Data.mat']])