function [Efunvals,Jvals] = evalCostFun_TimeDependentKoopEfun(Efun,y,u,Q,R,REF,t)

Efunvals = zeros(length(y),1);
Jvals = zeros(length(y),1);
for k=1:length(y)
    Efunvals(k) = Efun(y(k,:),t(k));
    Jvals(k) = Q*(Efunvals(k)-REF)^2 + u(k,:)*R*u(k,:)';
end


