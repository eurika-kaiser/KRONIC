function [Hvals1,Jvals1,uvals1] = evalCostFun_AsymPotentialWell_PhaseController(y1,B,R,Q,H,gain,a,REF)

Nc = size(B,2); Nt = size(y1,1);
uvals1 = zeros(Nc,Nt);
Hvals1 = zeros(Nt,1);
Jvals1 = zeros(Nt,1);

for i = 1:Nt
    Hvals1(i) = H(y1(i,:)); 
    % PHASE 1: control on in left well to drive to homoclinic orbit
    if  Hvals1(i) ~= H([a,0]) && y1(i,1) <= a
        uvals1(:,i) = -gain(y1(i,:))*(Hvals1(i)-H([a,0]));
%         HREF = H([a,0]);
        
        % PHASE 2: control off on homoclinic orbit while x <= a
    elseif  Hvals1(i) == H([a,0]) && y1(i,1) <= a
        uvals1(:,i) = zeros(Nc,1);
%         HREF = H([a,0]);
        
        % PHASE 3: control on in right well to drive to energy minimum
    elseif y1(i,1) > a
        uvals1(:,i) = -gain(y1(i,:))*(Hvals1(i)-H([1;0]));
    end
    Jvals1(i) = Q*(Hvals1(i)-REF)^2 + uvals1(:,i)'*R*uvals1(:,i);
end