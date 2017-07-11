function dx = AsymmetricDoubleWellFUN(t,x,gain,H,a,B)
Nc = size(B,2);

% PHASE 1: control on in left well to drive to homoclinic orbit
if  H(x) ~= H([a,0]) && x(1) <= a
    u = -gain(x)*(H(x)-H([a,0]));         
    disp('1')

% PHASE 2: control off on homoclinic orbit while x <= a
elseif  H(x) == H([a,0]) && x(1) <= a
    u = zeros(Nc,1);
    disp('2')
     
% PHASE 3: control on in right well to drive to energy minimum
elseif x(1) > a
    u = -gain(x)*(H(x)-H([1;0]));    
    disp('3')
end

dx = [x(2); -(x(1).^3 - x(1) - a.*x(1).^2 + a)]+B*u;
