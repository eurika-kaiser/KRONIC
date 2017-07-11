function [Hvals,Jvals] = evalCostFun_Hamiltonian(H,y,u,Q,R,REF)

Hvals = zeros(length(y),1);
Jvals = zeros(length(y),1);
for k=1:length(y)
    Hvals(k) = H(y(k,:));
    Jvals(k) = Q*(Hvals(k)-REF)^2 + u(:,k)'*R*u(:,k);
end


