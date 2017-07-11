clear all, close all, clc

path2figs = '../Figures/';
path2data = '../Data/';
ModelName = 'HamiltonianDuffing_';

% Parameters
Q  = 1;
R  = 1;
REF = 0;
a = 0;%sqrt(1/5);
V = 1;

% Hamiltonian
[X,Y] = meshgrid([-2.5:0.01:2.5], [-2.5:0.01:2.5]);
Hfield = (1/2)*Y.^2-(X.^2)/2 + (V/4)*X.^4;% 

% Initial condition 
y0ic    = [0 0];
xvec    = -1;
yvec    = -2.5:0.25:-0.25;
x       = [-4:0.01:4];
dt      = 0.001;
Ny      = length(yvec);
Nx      = length(xvec);
yIC     = zeros(2,Ny,Nx);
[x0,y0] = meshgrid(xvec+y0ic(1),yvec+y0ic(2));
yIC(1,:,:) = x0;
yIC(2,:,:) = y0;
yIC = [yIC, [1; 0.5], [1; 0.25]]; Ny = Ny+2;

%% UNFORCED
B       = [0; 0];
a       = 0;
V       = 1;
duration = 25;
tspan   =[0,duration]; 
L       = duration/dt;
yin     = yIC; 
yout    = zeros(2,Ny,Nx,L);
yout(:,:,:,1) = yin;
for step = 1:L-1
    time = step*dt;
    yout(:,:,:,step+1) = rk4singlestep(@(t,y,u)Duffing_ensemble(t,y,0,1,1,0,0,B,u),dt,time,yin,0);
    yin = yout(:,:,:,step);   
end


%% KRONIC
B = [0; 1];
ModelName1 = [ModelName, 'B01_'];
% B = [1; 0];
% ModelName = [ModelName, 'B10_'];
Bc = B;

B = Bc;
f = @(t,x,u)([x(2); x(1)-V*x(1)^3]+B*u);
H = @(x)((1/2)*x(2).^2-(x(1).^2)/2 + (V/4)*x(1).^4);
gradH = @(x)([-x(1) + V*x(1)^3; x(2)]);
gain  = @(x)(lqr(0,(gradH(x)'*B),Q,R));

% yIC = squeeze(yout(1:2,1,1,1:1000:8000));
a = 2.5;
yIC = [-2 0 2 -1 1 -2  0  2 -1  1 -a     -a -a  a     a a; ...
        a a a  a a -a -a -a -a -a -2  0.001  2 -2 0.001 2];
Nic     = size(yIC,2);    
yin = yIC;
uin = zeros(1,Nic);
for i = 1:Nic
    uin(i) = -gain(yin(:,i))*H(yin(:,i));
end
yout_ctrl           = zeros(2,Nic,L);
yout_ctrl(:,:,1)    = yin;
for step = 1:L-1
    time = step*dt;
    yout_ctrl(:,:,step+1) = rk4singlestep(@(t,y,u)Duffing_ensemble(t,y,0,1,1,0,0,B,u),dt,time,yin,uin);
    yin = yout_ctrl(:,:,step+1); 
    uin = zeros(1,Nic);
    for i = 1:Nic
        uin(i) = -gain(yin(:,i))*H(yin(:,i));
    end
end

% Inner circles
yIC = [squeeze(yout(1:2,end,1,1:3000:10000)),squeeze(yout(1:2,end-2,1,1:3000:10000))];
Ny = size(yIC,2);
yIC(1,1:Ny/2) = 1+0.1*(yIC(1,1:Ny/2)-1);
yIC(2,1:Ny/2) = 0.1*yIC(2,1:Ny/2);
yIC(1,Ny/2+1:end) = -1+0.1*(yIC(1,Ny/2+1:end)+1);
yIC(2,Ny/2+1:end) = 0.1*yIC(2,Ny/2+1:end);
Nic = size(yIC,2);
yin = yIC;
uin = zeros(1,Nic);
for i = 1:Nic
    uin(i) = -gain(yin(:,i))*H(yin(:,i));
end
yout_ctrl_in           = zeros(2,Nic,L);
yout_ctrl_in(:,:,1)    = yin;
for step = 1:L-1
    time = step*dt;
    yout_ctrl_in(:,:,step+1) = rk4singlestep(@(t,y,u)Duffing_ensemble(t,y,0,1,1,0,0,B,u),dt,time,yin,uin);
    yin = yout_ctrl_in(:,:,step+1); 
    uin = zeros(1,Nic);
    for i = 1:Nic
        uin(i) = -gain(yin(:,i))*H(yin(:,i));
    end
end



%% SAVE RESULTS
save([path2data,[ModelName1,'Ensemble.mat']])

