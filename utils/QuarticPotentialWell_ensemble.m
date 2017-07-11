function dy = QuarticPotentialWell3D(t,y)
% y is a two dimensional state-vector
dy = [  y(2,:,:);
      - y(1,:,:).^3];