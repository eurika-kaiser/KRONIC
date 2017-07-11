function [ynew,StartIDX,outside_vec,switched_domain_vec,Hvec] = applyPeriodicBC(y,xlim,ylim,Psi)

N = size(y,1);
ynew = zeros(size(y));
StartIDX = [1];
outside = 0;
outside_vec = zeros(size(y,1),1);
switched_domain_vec = zeros(size(y,1),1);
Hvec = zeros(size(y,1),1);
for i = 1:N
    ynew(i,:) = [mod(y(i,1),2),mod(y(i,2),2)];
    
        if (Psi(ynew(i,:))) > 0 && (ynew(i,2)>ylim(2)) % top right
            ynew(i,:) = ynew(i,:) + [-(xlim(2)-xlim(1))/2,-1];
        elseif (Psi(ynew(i,:))) < 0 && (ynew(i,2)>ylim(2)) % top left   
            ynew(i,:) = ynew(i,:) + [+(xlim(2)-xlim(1))/2,-1];
        end
    Hvec(i) = Psi(ynew(i,:));
    
    % Check whether outside or inside domain of interest
    if (y(i,1)>xlim(1)) && (y(i,1)<xlim(2)) && (y(i,2)>ylim(1)) && (y(i,2)<y(2))
        outside_new = 0;
    else
        outside_new = 1;
    end
    if outside ~= outside_new 
        switched_domain = 1;
    else
        switched_domain = 0;
    end
    outside = outside_new;
    outside_vec(i) = outside;
    switched_domain_vec(i) = switched_domain;
    
    while switched_domain == 1
        StartIDX = [StartIDX, i];
        switched_domain = 0;
    end
    
end
StartIDX = [StartIDX,N+1];


StartIDX = [1];
for i = 2:N
    dist = sum((ynew(i,:)-ynew(i-1,:)).^2);
    if dist>0.1
        StartIDX = [StartIDX,i];
    end
    
end
StartIDX = [StartIDX,N+1];
