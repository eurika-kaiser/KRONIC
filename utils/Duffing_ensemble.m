function dy = Duffing_ensemble(t,y,d,a,b,g,w,B,u)
% y is a two dimensional state-vector
dy = [  y(2,:,:) + +B(1)*u;
      - d*y(2,:,:) + a*y(1,:,:) - b*y(1,:,:).^3 + g*cos(w*t) + B(2)*u];
  