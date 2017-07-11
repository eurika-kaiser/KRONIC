function yout = buildGamma(yin,ydotin,nVars,polyorder,usesine)
% Copyright 2016, All Rights Reserved
% Code by Steven L. Brunton

n = size(yin,1);

ind = 1;
% poly order 0
% yout(:,ind) = ones(n,1);
% ind = ind+1;

% poly order 1
for i=1:nVars
    yout(:,ind) = ydotin(:,i);
    ind = ind+1;
end

if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            yout(:,ind) = ydotin(:,i).*yin(:,j) + yin(:,i).*ydotin(:,j);
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                yout(:,ind) = ydotin(:,i).*yin(:,j).*yin(:,k) + yin(:,i).*ydotin(:,j).*yin(:,k) + yin(:,i).*yin(:,j).*ydotin(:,k);
                ind = ind+1;
            end
        end
    end
end

if(polyorder>=4)
    % poly order 4
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
%                     yout(:,ind) = zeros(size(yin,1),1);
%                     ijkl = [i j k l];
%                     for k=1:4                         
%                         yout(:,ind) = yout(:,ind) + ydotin(:,ijkl(1)).*yin(:,ijkl(2)).*yin(:,ijkl(3)).*yin(:,ijkl(4));
%                         ijkl = [ijkl(2:end) ijkl(1)];
%                     end
                    yout(:,ind) = ydotin(:,i).*yin(:,j).*yin(:,k).*yin(:,l) + yin(:,i).*ydotin(:,j).*yin(:,k).*yin(:,l) +yin(:,i).*yin(:,j).*ydotin(:,k).*yin(:,l) +yin(:,i).*yin(:,j).*yin(:,k).*ydotin(:,l);
                    ind = ind+1;
                end
            end
        end
    end
end

if(polyorder>=5)
    % poly order 5
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                for l=k:nVars
                    for m=l:nVars
                        yout(:,ind) = ydotin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*yin(:,m) + yin(:,i).*ydotin(:,j).*yin(:,k).*yin(:,l).*yin(:,m) + yin(:,i).*yin(:,j).*ydotin(:,k).*yin(:,l).*yin(:,m) + yin(:,i).*yin(:,j).*yin(:,k).*ydotin(:,l).*yin(:,m) + yin(:,i).*yin(:,j).*yin(:,k).*yin(:,l).*ydotin(:,m);
                        ind = ind+1;
                    end
                end
            end
        end
    end
end

if(usesine)
    for k=1:10;
        yout = [yout sin(k*yin) cos(k*yin)];
    end
end