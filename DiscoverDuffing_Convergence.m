clear all, close all, clc
addpath('./utils');
path2figs = './../Figures/DUFFING/';
mkdir(path2figs)
ModelName = 'DiscoverDuffingCVG';

% Parameters
dt = 0.001;
tspan = dt:dt:10;
x0 = [0; -2.8];

% Function definitions
f = @(t,x)([x(2); x(1) - x(1)^3]);
H = @(x)( (1/2)*x(2).^2-(1/2)*x(1).^2 + (1/4)*x(1).^4 );
ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);

% Parameters for KRONIC
usesine = 0;
polyorder = 4;
nvar = 2;

%% Collect trajectory data
ModelName_tmp = [ModelName, '_trajectory_'];
[t,y] = ode45(f,tspan,x0,ode_options);

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

%% Compute derivatives
% using fourth order central difference
dy = zeros(length(y)-5,nvar);
for i=3:length(y)-3
    for k=1:nvar
        dy(i-2,k) = (1/(12*dt))*(-y(i+2,k)+8*y(i+1,k)-8*y(i-1,k)+y(i-2,k));
    end
end
x = y(3:end-3,1:nvar);
dx = dy;
tx = t(3:end-3);

%% Identify Koopman eigenfunction with lambda=0 (Hamiltonian) from full trajectory 
% Truth/reference
Hx = zeros(length(x),1);
for k=1:length(x)
    Hx(k) = H(x(k,:));
end

% Construct libraries
Theta = buildTheta(x,nvar,polyorder,usesine);
Gamma = buildGamma(x,dx,nvar,polyorder,usesine);

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
% xi0(abs(xi0)<1e-12) = 0;

D(IX(1))
xi1 = T(:,IX(1));%+Tls(:,IX(2));  % from least-squares fit
% xi1(abs(xi1)<1e-12) = 0; 

% Corrections for scaling
xi0 = 1/2*xi0/xi0(3); % energy
xi1 = 1/2*xi1/xi1(3); 

% Corrections for sign
if all((Theta)*(xi0)<0), xi0 = -xi0; end;
if all((Theta)*(xi1)<0), xi1 = -xi1; end;
     
% Print coefficients
poolDataLIST({'x','y'},xi0,nvar,polyorder,usesine);
poolDataLIST({'x','y'},xi1,nvar,polyorder,usesine);

% Plot evolution of eigenfunction = Hamiltonian
clear ph
figure; hold on, box on
ph(1) = semilogy(tx,((Theta)*(xi0))-Hx,'-b', 'LineWidth',2,'Color',[0,0,1]);
ph(2) = semilogy(tx,((Theta)*(xi1))-Hx,'-', 'LineWidth',2,'Color',[0,0.5,0]);
xlabel('t'), ylabel('E')
legend(ph,'SVD place', 'LS place')
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 260 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianReconstructionError','.eps']);

clear ph
figure; hold on, box on
ph(1) = semilogy(tx,Hx,'-k', 'LineWidth',18,'Color',[0.7,0.7,0.7]);
ph(2) = semilogy(tx,(Theta)*(xi0),'-b', 'LineWidth',8,'Color',[0,0,1]);
ph(3) = semilogy(tx,(Theta)*(xi1),'-', 'LineWidth',8,'Color',[0,0.5,0]);
xlabel('t'), ylabel('E')
legend(ph,'True place', 'SVD place', 'LS place')
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 260 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianReconstruction','.eps']);

% Plot error
clear ph
figure; hold on, box on
ph(1) = plot(tx,Gamma*xi0,'-k', 'LineWidth',2);%,'Color',[0,0,1]);
ph(2) = plot(tx,Gamma*xi1,'--r', 'LineWidth',2);%,'Color',[0,0.5,0]);
xlabel('t'), ylabel('Gamma xi')
xlim([min(t) max(t)])
ylim([min(Gamma*xi1)+0.2*min(Gamma*xi1) max(Gamma*xi1)+0.5*max(Gamma*xi1)])
legend(ph,'SVD p', 'LS  p')
set(gca,'xtick',[0:2:10])
set(gca,'FontSize',16)
set(gcf,'Position',[100 100 225 200])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'HamiltonianError','.eps']);

%% Error for increasing number of measurements
M = length(tx);
Theta0 = Theta;
Gamma0 = Gamma;

mset = [7,10,15,21,28,42,56,86,112,224,448,896, 1792, 3584, M];
Nm = length(mset);
HRerr_SVD = zeros(Nm,1);
HRerr_LSE = zeros(Nm,1);
HRerr_rel_SVD = zeros(Nm,1);
HRerr_rel_LSE = zeros(Nm,1);
Times_SVD = zeros(Nm,1);
Times_LSE = zeros(Nm,1);
HR_SVD = zeros(size(Theta0,1),Nm);
HR_LSE = zeros(size(Theta0,1),Nm);
NS_SVD = zeros(Nm,1);
NS_LSE = zeros(Nm,1);
Cost_SVD = zeros(Nm,1);
Cost_LSE = zeros(Nm,1);
Xc = zeros(size(y,1),size(y,2),Nm);
Hc = zeros(size(y,1),Nm);
uc = zeros(size(y,1),Nm);
xi_SVD = zeros(14,Nm);
xi_LSE = zeros(14,Nm);
%%
for iM = 1:length(mset)
    m = mset(iM);
    t0 = tic;
    
    % Data
    noise = 0*randn(m,nvar);
    ydat = y(1:m,:)+noise;
    
    tpre0 = tic;
    % Time derivative
    dy = zeros(length(ydat)-5,nvar);
    for i=3:length(ydat)-3
        for k=1:nvar
            dy(i-2,k) = (1/(12*dt))*(-ydat(i+2,k)+8*ydat(i+1,k)-8*ydat(i-1,k)+ydat(i-2,k));
        end
    end
    x = ydat(3:end-3,1:nvar);
    dx = dy;
    
    % Construct libraries
    Theta = buildTheta(x,nvar,polyorder,usesine);
    Gamma = buildGamma(x,dx,nvar,polyorder,usesine);
    tpre = toc(tpre0);
    
    tSVD = tic;
    % Compute SVD
    [U,S,V] = svd(0*Theta - Gamma,'econ');
    
    % Compute eigenfunction
    xi0 = V(:,end);             % from SVD
    xi0(abs(xi0)<1e-8) = 0; xi0 = xi0./norm(xi0);
    
    % Time to estimate
    Times_SVD(iM) = toc(tSVD)+tpre;
    
    tLSE = tic;
    % Least-squares Koopman
    K = pinv(Theta)*Gamma;
%     K(abs(K)<1e-8) = 0;
    [T,D] = eig(K);
    D = diag(D);
    [~,IX] = sort(abs(D),'ascend');
    
    if D(IX(1))==conj(D(IX(2)))
        xi1 = T(:,IX(1))+T(:,IX(2));
    else
        D(IX(1));
        xi1 = T(:,IX(1)); 
    end
    xi1(abs(xi1)<1e-8) = 0;
    xi1 = xi1./norm(xi1);
    
    Times_LSE(iM) = toc(tLSE)+tpre;
    
    xi_SVD(:,iM) = xi0;
    xi_LSE(:,iM) = xi1;
    
    % Print coefficients
%     poolDataLIST({'x','y'},xi0,nvar,polyorder,usesine);
%     poolDataLIST({'x','y'},xi1,nvar,polyorder,usesine);

    % Scaling corrections: Not necessary for control, just for estimating error 
    % Corrections for sign
    if all((Theta0)*(xi0)<0), xi0 = -xi0; end;
    if all((Theta0)*(xi1)<0), xi1 = -xi1; end;
    % Corrections for scaling
    xi0 = 1/2*xi0/xi0(3); % energy
    
    % Reconstruction Error
    HRerr_SVD(iM) = max(abs( (Theta0*(xi0))-Hx ));
    HRerr_rel_SVD(iM) = max(abs( ((Theta0*(xi0))-Hx)./Hx));
    HRerr_LSE(iM) = max(abs( (Theta0*(xi1))-Hx ));
    HRerr_rel_LSE(iM) = max(abs( ((Theta0*(xi1))-Hx)./Hx));
    
    % Null solution error
    NS_SVD(iM) = max(abs(Gamma0*xi0));
    NS_LSE(iM) = max(abs(Gamma0*xi1));
    
    % Eigenfunction time series
    HR_SVD(:,iM) = (Theta0)*(xi0);
    HR_LSE(:,iM) = (Theta0)*(xi1);
    
    % Control with KRONIC
    Q  = 1;
    R  = 1;
    REF = 0;
    A = 0;
    B = [0; 1];
    Aphi = 0; % Koopman eigenvalue = 0
    Bphi = @(x) [ [(buildThetaGradient(x',1,nvar,polyorder,usesine)*xi0); (buildThetaGradient(x',2,nvar,polyorder,usesine)*xi0)]' * B];
    gain = @(x)(lqr(A,Bphi(x),Q,R));
    Hfun = @(x) [buildTheta(x',nvar,polyorder,usesine)*xi0];
    fun = @(t,x,u)([x(2); x(1) - x(1)^3] + B*u);
    ode_options = odeset('RelTol',1e-10, 'AbsTol',1e-11);
    [~,y1] = ode45(@(t,x)fun(t,x,-gain(x)*(Hfun(x)-REF)),tspan,x0); %ode_options
    
    for i = 1:size(y1,1) % Calculate eigenfunction value and applied control
        Hc(i,iM) = H(y1(i,:)); 
        uc(i,iM) = -gain(y1(i,:)')*(Hfun(y1(i,:)')-REF);
    end
    J = cumsum(Q*(Hc(:,iM)-REF).^2 + R*uc(:,iM).^2);
    Cost_SVD(iM) = J(end);
    
    Xc(:,:,iM) = y1;
    
    
    % Control using LSE solution
%     Bphi = @(x) [ [(buildThetaGradient(x',1,nvar,polyorder,usesine)*xi1); (buildThetaGradient(x',2,nvar,polyorder,usesine)*xi1)]' * B];
%     gain = @(x)(lqr(A,Bphi(x),Q,R));
%     Hfun = @(x) [buildTheta(x',nvar,polyorder,usesine)*xi1];
%     fun = @(t,x,u)([x(2); x(1) - x(1)^3] + B*u);
%     [~,y1] = ode45(@(t,x)fun(t,x,-gain(x)*(Hfun(x)-REF)),tspan,x0); %ode_options
%     
%     for i = 1:size(y1,1)
%         Hc(i,iM) = H(y1(i,:)); 
%         uc(i,iM) = -gain(y1(i,:)')*(Hfun(y1(i,:)')-REF);
%     end
%     J = cumsum(Q*(Hc(:,iM)-REF).^2 + R*uc(:,iM).^2);
%     Cost_LSE(iM) = J(end);
    
    
    disp(['FINI: ',num2str(iM), ' of ',num2str(Nm)])
    toc(t0) % one run
end

%% Kink in curve
mset(13)
xi_SVD(:,13)
NS_SVD(13)


%% Show results
clear ph
figure; 
semilogy(mset,HRerr_SVD,'-k', 'LineWidth',2);hold on, box on, axis tight
% semilogy(mset,HRerr_LSE,'--r', 'LineWidth',2);
% semilogy(mset(13),HRerr_SVD(13),'or')
% xlim([0,10^4])
xlabel('Points'), ylabel('RecErr')
set(gca,'FontSize',16,'ytick',[10^0,10^1,10^2,10^3])%,'xtick',[1000,5000,10000])
set(gcf,'Position',[100 100 225 100])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'DataPoints_Error','.eps']);

clear ph
figure; 
semilogy(mset,NS_SVD,'-k', 'LineWidth',2);hold on, box on, axis tight,
% semilogy(mset(13),NS_SVD(13),'or')
%ph(2) = 
% semilogy(mset,NS_LSE,'--r', 'LineWidth',2);
xlim([0,10^4])
xlabel('Points'), ylabel('NullSolErr')
set(gca,'FontSize',16,'ytick',[10^-9,10^-6,10^0],'xtick',[1000,5000,10000])
set(gcf,'Position',[100 100 240 105])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'DataPoints_NullError','.eps']);

clear ph
figure; 
ph(1) = plot(mset,Times_SVD,'-k', 'LineWidth',2);hold on, box on,axis tight
% ph(2) = plot(7:M,NS_LSE(7:M),'--r', 'LineWidth',2);
xlim([0,10^4])
xlabel('Points'), ylabel('Time')
set(gca,'FontSize',16,'ytick',[2*1e-3, 6*1e-3],'xtick',[1000,5000,10000])
set(gcf,'Position',[100 100 225 115])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'DataPoints_Times','.eps']);

clear ph
figure; 
semilogy(mset,Cost_SVD,'-k', 'LineWidth',2);hold on, box on, axis tight
xlabel('Points'), ylabel('Cost')
xlim([0,10^4])
set(gca,'FontSize',16,'ytick',[10^4,10^6],'xtick',[1000,5000,10^4])
set(gcf,'Position',[100 100 305 120])
set(gcf,'PaperPositionMode','auto')
print('-depsc2', '-loose', [path2figs,ModelName_tmp,'DataPoints_ControlCost','.eps']);

%% Show controlled trajectories
figure,hold on,box on
cmap = jet(Nm);
x1 = -sqrt(2):0.01:sqrt(2);
plot(x1,x1.*sqrt(1-0.5.*(x1).^2),'--k')
plot(x1,-x1.*sqrt(1-0.5.*(x1).^2),'--k')
iM = Nm; plot(Xc(:,1,iM),Xc(:,2,iM),'-','Color',cmap(iM,:),'LineWidth',4);
iM = Nm-1; plot(Xc(:,1,iM),Xc(:,2,iM),'-','Color',cmap(iM,:),'LineWidth',3);
iM = Nm-2; plot(Xc(:,1,iM),Xc(:,2,iM),'-','Color',cmap(iM,:),'LineWidth',2);
for iM = Nm:-1:1
    plot(Xc(:,1,iM),Xc(:,2,iM),'-','Color',cmap(iM,:),'LineWidth',3-(Nm-iM)*0.2);
end
set(gcf,'Position',[100 100 305 120])
set(gcf,'PaperPositionMode','auto')


