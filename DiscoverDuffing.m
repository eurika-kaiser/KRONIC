clear all, close all, clc
addpath('./utils');
% addpath('../../../Sources/Code/psv-master')
path2figs = './../Figures/DUFFING/';
mkdir(path2figs)
ModelName = 'DiscoverDuffing';

% Parameters
dt = 0.001;
tspan = dt:dt:10;
x0 = [0; -2.8];

f = @(t,x)([x(2); x(1) - x(1)^3]);
H = @(x)( (1/2)*x(2).^2-(1/2)*x(1).^2 + (1/4)*x(1).^4 );
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);

usesine = 0;
polyorder = 4;
nvar = 2;

%% Trajectory data
ModelName_tmp = [ModelName, '_trajectory_'];

[t,y] = ode45(f,tspan,x0,ode_options);

Hy = zeros(length(y),1);
dy = zeros(size(y));
for k=1:length(y)
    dy(k,:) = f(0,y(k,:))';
    Hy(k) = H(y(k,:));
end

figure; hold on, box on
plot(t,y(:,1),'-','Color',[0,0,0.7],'LineWidth',2)
plot(t,y(:,2),'-','Color',[0,0.7,0],'LineWidth',2)
legend('x1','x2')
xlabel('t'), ylabel('xi')
set(gca,'xtick',[0:2:10])
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'Hamiltonian_Trajectory','.eps']);

% Construct libraries
Theta = buildTheta(y,nvar,polyorder,usesine);
Gamma = buildGamma(y,dy,nvar,polyorder,usesine);

% Compute SVD
[U,S,V] = svd(0*Theta - Gamma,'econ');

% Least-squares Koopman
K = pinv(Theta)*Gamma;
K(abs(K)<1e-12) = 0;
[T,D] = eig(K);
D = diag(D);
[~,IX] = sort(abs(D),'ascend');

% Compute eigenfunction
xi0 = V(:,end);             % from SVD
xi0(abs(xi0)<1e-12) = 0;

D(IX(1))
xi1 = T(:,IX(1));%+Tls(:,IX(2));  % from least-squares fit
xi1(abs(xi1)<1e-12) = 0; 

% Print coefficients
poolDataLIST({'x','y'},xi0,nvar,polyorder,usesine);
poolDataLIST({'x','y'},xi1,nvar,polyorder,usesine);

% Plot evolution of eigenfunction = Hamiltonian
if length(Hy)~=length(t)
    t = 1:length(Hy);
end

figure; hold on, box on
ph(1) = plot(t,Hy./norm(Hy),'-k', 'LineWidth',18,'Color',[0.7,0.7,0.7]);
ph(2) = plot(t,(Theta)*(xi0)./norm((Theta)*(xi0)),'-b', 'LineWidth',8,'Color',[0,0,1]);
ph(3) = plot(t,-(Theta)*(xi1)./norm((Theta)*(xi1)),'--', 'LineWidth',8,'Color',[0,0.5,0]);
xlabel('t'), ylabel('E')
ylim([0.01-1e-4 0.01+2e-4]), xlim([min(t) max(t)])
set(gca,'xtick',[0:2:10])
legend(ph,'True place', 'SVD place', 'LS place')
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 260 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianReconstruction','.eps']);

% Plot error 
dstep = 5;
clear ph
figure; hold on, box on
tmp = Gamma*xi0;
ph(1) = plot(t(1:dstep:end),tmp(1:dstep:end),'-k', 'LineWidth',2);%,'Color',[0,0,1]);
tmp = -Gamma*xi1;
ph(2) = plot(t(1:dstep:end),tmp(1:dstep:end),'-r', 'LineWidth',2);%,'Color',[0,0.5,0]);
xlabel('t'), ylabel('Gamma xi')
xlim([min(t) max(t)])
ylim([min(Gamma*xi1)+0.2*min(Gamma*xi1) max(Gamma*xi1)+0.5*max(Gamma*xi1)])
legend(ph,'SVD p', 'LS  p')
set(gca,'xtick',[0:2:10])
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianError','.eps']);

err0 = (Hy./norm(Hy)-(Theta)*(xi0)./norm((Theta)*(xi0)))./(Hy./norm(Hy));
err1 = (Hy./norm(Hy)-(Theta)*(xi1)./norm((Theta)*(xi1)))./(Hy./norm(Hy));
figure, plot(err0,'-r'),hold on, plot(err1,'--b')
axis tight

return
%% Ensemble data
close all
ModelName_tmp = [ModelName, '_ensemble_'];

% Collect data
clear y dy Hy
[Xgrid,Vgrid] = meshgrid([-2.45:.05:2.5],[-2.45:.05:2.5]);
Nx = size(Xgrid,1), Ny = size(Vgrid,2);
y(:,1) = reshape(Xgrid,[Nx*Ny 1]);
y(:,2) = reshape(Vgrid,[Nx*Ny 1]);
Hy = zeros(length(y),1);
Hdy = zeros(length(y),1);
dy = zeros(size(y));
for k=1:length(y)
    dy(k,:) = f(0,y(k,:))';
    Hy(k) = H(y(k,:));
    Hdy(k) = H(dy(k,:));
end

cmap = gray(2*11);
cmap = cmap(1:2:end,:);
xH = [-2.4:0.01:2.4];
[X,V] = meshgrid(xH,xH);
Hfield = (1/2)*V.^2-(X.^2)/2 + (1/4)*X.^4; 

figure; box on, hold on
sh = surf(xH,xH,zeros(size(Hfield)));shading interp, view(2)
sh.FaceColor = [0 0 0];
contourf(xH,xH,log(Hfield+(0.2500001)),[log([-0.2,-0.1,0,0.25:0.5:4]+(0.2500001))],'LineColor','none'), colormap(cmap)
hold on
sh = scatter(Xgrid(:),Vgrid(:),'.b');
sh.SizeData = 5;
x1 = -sqrt(2):0.01:sqrt(2);
plot(x1,x1.*sqrt(1-0.5.*(x1).^2),'--y')
plot(x1,-x1.*sqrt(1-0.5.*(x1).^2),'--y')
xlabel('x1'), ylabel('x2')
drawnow
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 220 200])
set(gcf,'PaperPositionMode','auto')
% print('-painters','-depsc2', '-loose', [path2figs,ModelName_tmp,'Hamiltonian_SamplingPoints','.eps']);
print('-opengl','-depsc2', '-loose', [path2figs,ModelName_tmp,'Hamiltonian_SamplingPoints','.eps']);

% Construct libraries
Theta = buildTheta(y,nvar,polyorder,usesine);
Gamma = buildGamma(y,dy,nvar,polyorder,usesine);

% Compute SVD
[U,S,V] = svd(0*Theta - Gamma,'econ');

% Least-squares Koopman
K = pinv(Theta)*Gamma;
K(abs(K)<1e-12) = 0;
[T,D] = eig(K);
D = diag(D);
[~,IX] = sort(abs(D),'ascend');

% Compute eigenfunction
xi0 = V(:,end);             % from SVD
xi0(abs(xi0)<1e-12) = 0;

D(IX(1))
xi1 = T(:,IX(1));  % from least-squares fit
xi1(abs(xi1)<1e-12) = 0; 

% Print coefficients
poolDataLIST({'x','y'},xi0,nvar,polyorder,usesine);
poolDataLIST({'x','y'},xi1,nvar,polyorder,usesine);

% Plot evolution of eigenfunction = Hamiltonian
if length(Hy)~=length(t)
    t = 1:length(Hy);
end

%% Show results for ensemble of points
close all
figure; hold on, box on
ph(1) = plot(t,Hy./norm(Hy),'-k', 'LineWidth',5,'Color',[0.7,0.7,0.7]);
ph(2) = plot(t,(Theta)*(xi0)./norm((Theta)*(xi0)),'-k', 'LineWidth',2);%,'Color',[0,0,1]);
ph(3) = plot(t,-(Theta)*(xi1)./norm((Theta)*(xi1)),'-r', 'LineWidth',2);%,'Color',[0,0.5,0]);
xlabel('sampling point'), ylabel('E')
ylim([-0.005 0.05]),
xlim([-0 10])
legend(ph,'True place', 'SVD place', 'LS place','location','north')
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 240 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianReconstruction','.eps']);

% Plot error
dstep = 5;
clear ph
figure; hold on, box on
tmp = Gamma*xi1;
ph(2) = plot([1:dstep:Nx*Ny],tmp(1:dstep:end),'-r', 'LineWidth',2);%,'Color',[0,0.5,0]);
tmp = Gamma*xi0;
ph(1) = plot([1:dstep:Nx*Ny],tmp(1:dstep:end),'-k', 'LineWidth',2);%,'Color',[0,0,1]);
xlabel('sampling point'), ylabel('Gamma xi')
axis tight
xlim([0,Nx*Ny+10])
set(gca, 'xtick',[5000 10000],'xticklabel',{'5k', '10k'})
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianError','.eps']);
