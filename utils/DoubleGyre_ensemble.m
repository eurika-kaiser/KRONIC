function dy = DoubleGyre_ensemble(t,y,A,omega,epsilon,B,u)
% y is a two dimensional state-vector
f = epsilon*sin(omega*t)*y(1,:,:).^2 + (1-2*epsilon*sin(omega*t)).*y(1,:,:);
dfdx = 2*epsilon*sin(omega*t)*y(1,:,:) + (1-2*epsilon*sin(omega*t));

if size(B,2) < 2
dy = [-A*pi*sin(pi*f).*cos(pi*y(2,:,:))       + B(1,:)*u(:,:,:); 
      A*pi*cos(pi*f).*sin(pi*y(2,:,:)).*dfdx  + B(2,:)*u(:,:,:)];
elseif size(B,2) >= 2  
dy = [-A*pi*sin(pi*f).*cos(pi*y(2,:,:)); 
      A*pi*cos(pi*f).*sin(pi*y(2,:,:)).*dfdx];
  
db = zeros(size(dy));  
for i = 1:size(B,2)
    db = db + [B(1,i)*u(i,:,:); B(2,i)*u(i,:,:)];
end

dy = dy+db;

end

          
  