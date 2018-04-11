function J = ObjectiveFCN_EDMDc(uopt,x,N,Nu,xref,u0,p,Q,R,Ru)
%% Cost function of nonlinear MPC
%
% Inputs:
%   u:      optimization variable, from time k to time k+N-1
%   x:      current state at time k
%   Ts:     controller sample time
%   N:      prediction horizon
%   xref:   state references, constant from time k+1 to k+N
%   u0:     previous controller output at time k-1
%
% Output:
%   J:      objective function cost
%

%% Nonlinear MPC design parameters
u = uopt;

%% Integrate system
[xk,~] = lsim(p.sys,[u' 0],[0:N].*p.sys.Ts,poolData(x',2,p.polyorder,p.usesine)');
xk = xk';
xk = xk(1:p.nvar,:);

%% Cost Calculation
% Set initial plant states, controller output and cost.
uk = u(1);
J = 0;
% Loop through each prediction step.
for ct=1:N
    % Obtain plant state at next prediction step.
    xk1 = xk(:,ct);

    % accumulate state tracking cost from x(k+1) to x(k+N).
    J = J + (p.Hfun(xk1)-p.Hfun(xref))'*Q*(p.Hfun(xk1)-p.Hfun(xref));
    % accumulate MV rate of change cost from u(k) to u(k+N-1).
    if ct==1
        J = J + (uk-u0)'*R*(uk-u0) + uk'*Ru*uk;
    else
        J = J + (uk-u(ct-1))'*R*(uk-u(ct-1)) + uk'*Ru*uk;
    end
    % Update uk for the next prediction step.
    if ct<N
        uk = u(ct+1);
    end
end


