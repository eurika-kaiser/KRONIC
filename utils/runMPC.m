function [xHistory, uHistory, tHistory, rHistory] = runMPC(syshandle, Duration,Ts,N,Nu,x0, ObjectiveFCN, ConstraintFCN, Q, R, Ru, options, xref, pest)
% Prepare variables
Nvar = length(xref);
Nt = (Duration/Ts)+1;
uopt0    = 0;
xhat     = x0;
uopt     = uopt0.*ones(Nu,1);
xHistory = zeros(Nvar,Nt); xHistory(:,1) = xhat;
uHistory = zeros(1,Nt);    uHistory(1)   = uopt(1);
tHistory = zeros(1,Nt);    tHistory(1)   = 0;
rHistory = zeros(Nvar,Nt);

% Start simulation
fprintf('Simulation started.  It might take a while...\n')
tic
for ct = 1:Nt-1
    
    % NMPC with full-state feedback
    COSTFUN = @(u) ObjectiveFCN(u,xhat,N,Nu,xref,uHistory(:,ct),pest,Q,R,Ru);
%     CONSFUN = @(u) ConstraintFCN(u,uHistory(:,ct),xhat,N,LBo,UBo,LBdu,UBdu,pest);
    uopt = fmincon(COSTFUN,uopt,[],[],[],[],[],[],[],options);
%     uopt = fmincon(COSTFUN,uopt,[],[],[],[],LB,UB,CONSFUN,option);
    
    % Integrate system
    xhat = rk4u(syshandle,xhat,uopt(1),Ts/10,10,[],0); %10, 2
    xHistory(:,ct+1) = xhat;
    uHistory(:,ct+1) = uopt(1);
    tHistory(:,ct+1) = ct*Ts;
    rHistory(:,ct+1) = xref;
    
    if mod(ct,10) == 0
        disp(['PROGRESS: ',num2str(100*ct/(Duration/Ts)),'%'])
    end
end