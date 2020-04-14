clear all, close all, clc
addpath('./utils');
path2figs = './../Figures/DUFFING/'; mkdir(path2figs)
ModelName = 'DiscoverDuffing_KRONICvsEDMDc';

% Parameters
dt = 0.01;
tspan = dt:dt:10;
x0 = [0; -2.8];

% Function definitions
f = @(t,x)([x(2); x(1) - x(1)^3]);
H = @(x)( (1/2)*x(2).^2-(1/2)*x(1).^2 + (1/4)*x(1).^4 );
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
nvar = 2;

%% Collect trajectory data
% Unforced
[t,y] = ode45(f,tspan,x0,ode_options);

% Forced
B = [0; 1];
forcing = @(x,t) [(2*(sin(2*pi*.1*t).*sin(2*pi*1*t)))];
f = @(t,x) [x(2); x(1) - x(1)^3]  + B*forcing(x,t); 
[tF,yF] = ode45(f,tspan,x0,ode_options);

uF = forcing(yF,tspan)';
U = uF';

figure; hold on, box on
plot3(y(:,1),y(:,2),zeros(size(y(:,2))),'-','Color',[0,0,0.7],'LineWidth',2)
color_line3(yF(:,1),yF(:,2),zeros(size(yF(:,2))),uF,'LineWidth',2);%,'--','Color','r','LineWidth',2)
xlabel('x1'), ylabel('x2')

figure,plot(uF)

%% EDMDc
% Parameters
pEDMDc.usesine = 0;
pEDMDc.polyorder = 4;
Y = poolData(yF,2,pEDMDc.polyorder,pEDMDc.usesine)';
Nstates = size(Y,1);

G = [Y(:,1:end-1);U(:,1:end-1)];
[Ug,Sg,Vg] = svd(G,'econ');
AB = Y(:,2:end)*Vg*Sg^(-1)*Ug';
Am = AB(1:Nstates,1:Nstates);
Bm = AB(1:Nstates,end);
Cm = eye(Nstates,Nstates);
Dm = zeros(Nstates,1);
sysmodel_EDMDc = ss(Am,Bm,Cm,Dm,dt);

% yout = poolDataLIST(yin,ahat,nVars,polyorder,usesine)
yEDMDc = lsim(sysmodel_EDMDc,uF,tspan,poolData(x0',2,pEDMDc.polyorder,pEDMDc.usesine)');
yEDMDc = yEDMDc(:,1:2);

figure; hold on, box on
plot3(y(:,1),y(:,2),zeros(size(y(:,2))),'-','Color',[0,0,0.7],'LineWidth',2)
color_line3(yF(:,1),yF(:,2),zeros(size(yF(:,2))),uF,'LineWidth',2);%,'--','Color','r','LineWidth',2)
plot3(yEDMDc(:,1),yEDMDc(:,2),zeros(size(yEDMDc(:,2))),'--','Color','r','LineWidth',2)
xlabel('x1'), ylabel('x2')
legend('Unforced','Forced','EDMDc')
%% KRONIC
% Parameters
pKRONIC.usesine = 0;
pKRONIC.polyorder = 4;

% Compute derivatives using fourth order central difference
dy = zeros(length(y)-5,nvar);
for i=3:length(y)-3
    for k=1:nvar
        dy(i-2,k) = (1/(12*dt))*(-y(i+2,k)+8*y(i+1,k)-8*y(i-1,k)+y(i-2,k));
    end
end
x = y(3:end-3,1:nvar);
dx = dy;
tx = t(3:end-3);

% Construct libraries
Theta = buildTheta(x,nvar,pKRONIC.polyorder,pKRONIC.usesine);
Gamma = buildGamma(x,dx,nvar,pKRONIC.polyorder,pKRONIC.usesine);

% Compute SVD
[U,S,V] = svd(0*Theta - Gamma,'econ');
KRONIC.phi = V(:,end);

% Remove values below small threshold
KRONIC.phi(abs(KRONIC.phi)<1e-6) = 0;

% Print coefficients
poolDataLIST({'x','y'},KRONIC.phi,nvar,pKRONIC.polyorder,pKRONIC.usesine);

%% Control parameters & functions 
x_REF = [0; 0];
H_REF = H(x_REF);

f = @(t,x,u,p)([x(2); x(1)-x(1)^3] + B*u);
Hc = @(x)((1/2)*x(2).^2-(x(1).^2)/2 + (1/4)*x(1).^4 - H_REF);
gradH = @(x)([-x(1) + x(1)^3; x(2)]);

% Parameters MPC
options = optimoptions('fmincon','Algorithm','sqp','Display','none', ...
    'MaxIterations',100);
Duration = 10;                  % Run for 'Duration' time units
N  = 5;                         % Prediction horizon (number of iterations)
Nu  = N;                        % Control horizon (number of iterations)
Q = 5*eye(2);                   % x weights
Ru = 1;                         % u weights
R = 0;                          % du weights
QH = 5;                         % phi(=H)  weights

%% EDMDc 
% MPC
p_EDMDc.sys = sysmodel_EDMDc;
p_EDMDc.nvar = 2;
p_EDMDc.polyorder = pEDMDc.polyorder;
p_EDMDc.usesine = pEDMDc.usesine; 
p_EDMDc.Hfun = H;
tic
[xHistory_EDMDc_MPC, uHistory_EDMDc_MPC, tHistory_EDMDc_MPC, rHistory_EDMDc_MPC] = runMPC(f,Duration,dt,N,Nu,x0,@ObjectiveFCN_EDMDc,[],QH,R,Ru,options,x_REF,p_EDMDc);
toc

% EDMDc // LQR fails
% Qlarge = zeros(Nstates,Nstates); Qlarge(1:nvar,1:nvar) = Q;
% K = lqr(sysmodel_EDMDc.A,sysmodel_EDMDc.B,Qlarge,Ru);
% [~,yEDMDc] = ode45(@(t,x)f(t,x,-K*poolData(x',2,polyorder,usesine)'),tspan,x0);

%% KRONIC (analytic)

% MPC
Hsys = @(t,x,u,p) [gradH(x)'*B*u];
% p_KRONIC.sys = Hsys;
p_analyticKRONIC.A = 0;
p_analyticKRONIC.B = @(x) [gradH(x)'*B];
p_analyticKRONIC.H = H;
p_analyticKRONIC.Ts = dt;
tic
[xHistory_analyticKRONIC_MPC, uHistory_analyticKRONIC_MPC, tHistory_analyticKRONIC_MPC, rHistory_analyticKRONIC_MPC] = runMPC(f,Duration,dt,N,Nu,x0,@ObjectiveFCN_KRONIC,[],QH,R,Ru,options,x_REF,p_analyticKRONIC);
toc

% LQR
tic
gain = @(x)(lqr(0,(gradH(x)'*B),QH,Ru));
[tHistory_analyticKRONIC_LQR,xHistory_analyticKRONIC_LQR] = ode45(@(t,x)f(t,x,-gain(x)*Hc(x)),tspan,x0);
xHistory_analyticKRONIC_LQR = xHistory_analyticKRONIC_LQR';
uHistory_analyticKRONIC_LQR = zeros(size(tHistory_analyticKRONIC_LQR));
for i = 1:length(xHistory_analyticKRONIC_LQR)
    uHistory_analyticKRONIC_LQR(i) = -gain(xHistory_analyticKRONIC_LQR(:,i))*Hc(xHistory_analyticKRONIC_LQR(:,i));
end
toc

%% KRONIC (data)

% MPC
pKRONIC.phi = KRONIC.phi;
pKRONIC.A = 0;
pKRONIC.B = @(x) [(buildThetaGradient(x',1,nvar,pKRONIC.polyorder,pKRONIC.usesine)*pKRONIC.phi);(buildThetaGradient(x',2,nvar,pKRONIC.polyorder,pKRONIC.usesine)*pKRONIC.phi)]'*B;
pKRONIC.H = @(x) buildTheta(x',nvar,pKRONIC.polyorder,pKRONIC.usesine)*pKRONIC.phi;
pKRONIC.Ts = dt;
tic
[xHistory_KRONIC_MPC, uHistory_KRONIC_MPC, tHistory_KRONIC_MPC, rHistory_KRONIC_MPC] = runMPC(f,Duration,dt,N,Nu,x0,@ObjectiveFCN_KRONIC,[],QH,R,Ru,options,x_REF,pKRONIC);
toc


% LQR
tic
pKRONIC.Hc = @(x) [buildTheta(x',nvar,pKRONIC.polyorder,pKRONIC.usesine)*KRONIC.phi - pKRONIC.H(x_REF)];   
gain = @(x)(lqr(0,pKRONIC.B(x),QH,Ru));
[tHistory_KRONIC_LQR,xHistory_KRONIC_LQR] = ode45(@(t,x)f(t,x,-gain(x)*pKRONIC.Hc(x)),tspan,x0);
xHistory_KRONIC_LQR = xHistory_KRONIC_LQR';
uHistory_KRONIC_LQR = zeros(size(tHistory_KRONIC_LQR));
for i = 1:length(xHistory_KRONIC_LQR)
    uHistory_KRONIC_LQR(i) = -gain(xHistory_KRONIC_LQR(:,i))*pKRONIC.Hc(xHistory_KRONIC_LQR(:,i));
end
toc


%%
figure; hold on, box on
plot3(y(:,1),y(:,2),zeros(size(y(:,2))),'-','Color',[0,0,0.7],'LineWidth',2)
% color_line3(yF(:,1),yF(:,2),zeros(size(yF(:,2))),uF,'LineWidth',2);%,'--','Color','r','LineWidth',2)
% plot3(yEDMDc(:,1),yEDMDc(:,2),zeros(size(yEDMDc(:,2))),'--','Color','r','LineWidth',2)
% color_line3(xHistory(1,:),xHistory(2,:),zeros(size(xHistory(1,:))),uHistory','LineWidth',2);%,'--','Color','r','LineWidth',2)
plot3(xHistory_EDMDc_MPC(1,:),xHistory_EDMDc_MPC(2,:),zeros(size(xHistory_EDMDc_MPC(1,:))),'--','Color','g','LineWidth',2)
plot3(xHistory_analyticKRONIC_MPC(1,:),xHistory_analyticKRONIC_MPC(2,:),zeros(size(xHistory_analyticKRONIC_MPC(1,:))),'-','Color','b','LineWidth',2)
plot3(xHistory_KRONIC_MPC(1,:),xHistory_KRONIC_MPC(2,:),zeros(size(xHistory_KRONIC_MPC(1,:))),'--','Color','cyan','LineWidth',2)

plot3(xHistory_analyticKRONIC_LQR(1,:),xHistory_analyticKRONIC_LQR(2,:),zeros(size(xHistory_analyticKRONIC_LQR(1,:))),'-','Color','r','LineWidth',2)
plot3(xHistory_KRONIC_LQR(1,:),xHistory_KRONIC_LQR(2,:),zeros(size(xHistory_KRONIC_LQR(1,:))),'--','Color','m','LineWidth',2)

% color_line3(xHistory2(1,:),xHistory2(2,:),zeros(size(xHistory2(1,:))),uHistory2','LineWidth',2);%,'--','Color','r','LineWidth',2)
x1 = -sqrt(2):0.01:sqrt(2);
plot(x1,x1.*sqrt(1-0.5.*(x1).^2),'--k')
plot(x1,-x1.*sqrt(1-0.5.*(x1).^2),'--k')
xlabel('x'),ylabel('y')
legend('Unforced','EDMD-MPC','analytic KRONIC-MPC', 'KRONIC-MPC','analytic KRONIC-SDRE', 'KRONIC-SDRE','H=0')
axis tight
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 600 400])
set(gcf,'PaperPositionMode','auto')
%% Show Actuation
figure,plot(tHistory_EDMDc_MPC,uHistory_EDMDc_MPC,'-g','LineWidth',2), hold on, 
plot(tHistory_analyticKRONIC_MPC,uHistory_analyticKRONIC_MPC,'-b','LineWidth',2), 
plot(tHistory_KRONIC_MPC,uHistory_KRONIC_MPC,'--','Color','cyan','LineWidth',2),
plot(tHistory_analyticKRONIC_LQR,uHistory_analyticKRONIC_LQR,'-r','LineWidth',2),
plot(tHistory_KRONIC_LQR,uHistory_KRONIC_LQR,'--m','LineWidth',2),
xlabel('time'),ylabel('u')
legend('EDMD-MPC','analytic KRONIC-MPC', 'KRONIC-MPC','analytic KRONIC-SDRE', 'KRONIC-SDRE','H=0')
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 600 400])
set(gcf,'PaperPositionMode','auto')
%% Show eigenfunction time series (=energy/Hamiltonian)
Nt = min([length(tHistory_EDMDc_MPC),length(tHistory_KRONIC_MPC),length(tHistory_KRONIC_LQR)]);
H_EDMDc_MPC = zeros(Nt,1);
H_KRONIC_MPC = zeros(Nt,1);
H_KRONIC_LQR = zeros(Nt,1);
for i = 1:Nt
    H_EDMDc_MPC(i) = pKRONIC.H(xHistory_EDMDc_MPC(:,i));
    H_KRONIC_MPC(i) = pKRONIC.H(xHistory_KRONIC_MPC(:,i));
    H_KRONIC_LQR(i) = pKRONIC.H(xHistory_KRONIC_LQR(:,i));
end


figure,plot(H_EDMDc_MPC,'-k','LineWidth',2), hold on, plot(H_KRONIC_MPC,'--r','LineWidth',2), plot(H_KRONIC_LQR,':b','LineWidth',2)
xlabel('time'),ylabel('H')
legend('EDMD-MPC','KRONIC-MPC','KRONIC-LQR')
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 600 400])
set(gcf,'PaperPositionMode','auto')
%% Show cumulative cost in terms of state
LineWidth = 3;
[Jvals] = evalCostFun(xHistory_EDMDc_MPC(:,1:end-1),uHistory_EDMDc_MPC(:,1:end-1),Q,Ru,x_REF);
[Jvals2] = evalCostFun(xHistory_KRONIC_MPC(:,1:end-1),uHistory_KRONIC_MPC(:,1:end-1),Q,Ru,x_REF);
[JvalsKRONIC] = evalCostFun(xHistory_KRONIC_LQR,uHistory_KRONIC_LQR',Q,Ru,x_REF);
figure,box on
plot(tHistory_KRONIC_LQR,cumsum(Jvals),'--r','LineWidth',LineWidth), hold on, 
plot(tHistory_KRONIC_LQR,cumsum(Jvals2),':k','LineWidth',LineWidth), 
plot(tHistory_KRONIC_LQR,cumsum(JvalsKRONIC),'-b','LineWidth',LineWidth)
ylim([0 2.5e4])
ylabel('J'), xlabel('t')
set(gca,'FontSize',16)
% axis tight
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [path2figs,ModelName,'_PerformanceComparisonWithKRONIC_Jx','.eps']);

%% Show cumulative cost in terms of eigenfunction
clear ph
[Jvals_H] = evalCostFun(H_EDMDc_MPC',uHistory_EDMDc_MPC(:,1:end-1),QH,Ru,H(x_REF));
[Jvals2_H] = evalCostFun(H_KRONIC_MPC',uHistory_KRONIC_MPC(:,1:end-1),QH,Ru,H(x_REF));
[JvalsKRONIC_H] = evalCostFun(H_KRONIC_LQR',uHistory_KRONIC_LQR',QH,Ru,H(x_REF));
figure,box on
ph(1) = plot(tHistory_KRONIC_LQR,cumsum(Jvals_H),'--r','LineWidth',LineWidth); hold on, 
ph(2) = plot(tHistory_KRONIC_LQR,cumsum(Jvals2_H),':k','LineWidth',LineWidth);
ph(3) = plot(tHistory_KRONIC_LQR,cumsum(JvalsKRONIC_H),'-b','LineWidth',LineWidth);
ylim([0 2.5e4])
ylabel('J'), xlabel('t')
set(gca,'FontSize',16)
% axis tight
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-painters','-depsc2', '-loose', [path2figs,ModelName,'_PerformanceComparisonWithKRONIC_Jphi','.eps']);

axis off
ylim([5e4 1e5])
legend(ph,{'EDMDc-MPC','KRONIC-MPC', 'KRONIC-SDRE00'})
print('-painters','-depsc2', '-loose', [path2figs,ModelName,'_PerformanceComparisonWithKRONIC_legend','.eps']);

