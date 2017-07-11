function dy = SlowManifold_ensemble(t,y,mu,lambda)
% y is a three dimensional state-vector
dy = [  mu*y(1,:,:);
        lambda*( y(2,:,:) - y(1,:,:).^2)];