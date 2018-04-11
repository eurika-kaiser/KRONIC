function yout = buildThetaGradient(yin,iDenVar,nVars,polyorder,usesine)
n = size(yin,1);

ind = 1;
% poly order 0
% yout(:,ind) = ones(n,1);
% ind = ind+1;

% poly order 1
try
for i=1:nVars
    if i == iDenVar
        yout(:,ind) = ones(size(yin(:,i)));
    else
        yout(:,ind) = zeros(size(yin(:,i)));
    end
    ind = ind+1;
end
catch
    keyboard
end
if(polyorder>=2)
    % poly order 2
    for i=1:nVars
        for j=i:nVars
            if i == iDenVar && j == iDenVar
                yout(:,ind) = 2*yin(:,i);
            elseif i == iDenVar && j ~= iDenVar
                yout(:,ind) = yin(:,j);
            elseif i ~= iDenVar && j == iDenVar
                yout(:,ind) = yin(:,i);
            else
                yout(:,ind) = zeros(size(yin(:,i)));
            end
            ind = ind+1;
        end
    end
end

if(polyorder>=3)
    % poly order 3
    for i=1:nVars
        for j=i:nVars
            for k=j:nVars
                if i == iDenVar && j == iDenVar && k == iDenVar
                    yout(:,ind) = 3*yin(:,i).*yin(:,j);
                elseif i == iDenVar && j == iDenVar && k ~= iDenVar
                    yout(:,ind) = 2*yin(:,i).*yin(:,k);
                elseif i == iDenVar && j ~= iDenVar && k == iDenVar
                    yout(:,ind) = 2*yin(:,i).*yin(:,j);
                elseif i ~= iDenVar && j == iDenVar && k == iDenVar
                    yout(:,ind) = 2*yin(:,i).*yin(:,j);
                elseif i == iDenVar && j ~= iDenVar && k ~= iDenVar 
                    yout(:,ind) = yin(:,j).*yin(:,k);
                elseif i ~= iDenVar && j == iDenVar && k ~= iDenVar 
                    yout(:,ind) = yin(:,i).*yin(:,k);
                elseif i ~= iDenVar && j ~= iDenVar && k == iDenVar
                    yout(:,ind) = yin(:,i).*yin(:,j);
                else
                    yout(:,ind) = zeros(size(yin(:,k)));
                end
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
                    if i == iDenVar && j == iDenVar && k == iDenVar && l == iDenVar
                        yout(:,ind) = 4*yin(:,i).*yin(:,j).*yin(:,k);
                    elseif i == iDenVar && j == iDenVar && k == iDenVar && l ~= iDenVar
                        yout(:,ind) = 3*yin(:,i).*yin(:,j).*yin(:,l);
                    elseif i == iDenVar && j == iDenVar && k ~= iDenVar && l == iDenVar
                        yout(:,ind) = 3*yin(:,i).*yin(:,j).*yin(:,k);
                    elseif i == iDenVar && j ~= iDenVar && k == iDenVar && l == iDenVar  
                        yout(:,ind) = 3*yin(:,i).*yin(:,j).*yin(:,k);
                    elseif i ~= iDenVar && j == iDenVar && k == iDenVar && l == iDenVar 
                        yout(:,ind) = 3*yin(:,i).*yin(:,j).*yin(:,k);
                    elseif i == iDenVar && j == iDenVar && k ~= iDenVar && l ~= iDenVar   
                        yout(:,ind) = 2*yin(:,i).*yin(:,k).*yin(:,l);
                    elseif i == iDenVar && j ~= iDenVar && k == iDenVar && l ~= iDenVar 
                        yout(:,ind) = 2*yin(:,i).*yin(:,j).*yin(:,l);
                    elseif i ~= iDenVar && j == iDenVar && k == iDenVar && l ~= iDenVar 
                        yout(:,ind) = 2*yin(:,i).*yin(:,j).*yin(:,l);
                    elseif i == iDenVar && j ~= iDenVar && k ~= iDenVar && l == iDenVar 
                        yout(:,ind) = 2*yin(:,i).*yin(:,j).*yin(:,k);
                    elseif i ~= iDenVar && j == iDenVar && k ~= iDenVar && l == iDenVar 
                        yout(:,ind) = 2*yin(:,i).*yin(:,j).*yin(:,k);
                    elseif i ~= iDenVar && j ~= iDenVar && k == iDenVar && l == iDenVar  
                        yout(:,ind) = 2*yin(:,i).*yin(:,j).*yin(:,k);
                    elseif i == iDenVar && j ~= iDenVar && k ~= iDenVar && l ~= iDenVar
                        yout(:,ind) = yin(:,j).*yin(:,k).*yin(:,l);
                    elseif i ~= iDenVar && j == iDenVar && k ~= iDenVar && l ~= iDenVar
                        yout(:,ind) = yin(:,i).*yin(:,k).*yin(:,l);
                    elseif i ~= iDenVar && j ~= iDenVar && k == iDenVar && l ~= iDenVar 
                        yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,l);
                    elseif i ~= iDenVar && j ~= iDenVar && k ~= iDenVar && l == iDenVar 
                        yout(:,ind) = yin(:,i).*yin(:,j).*yin(:,k);
                    else
                        yout(:,ind) = zeros(size(yin(:,i)));
                    end
                    ind = ind+1;
                end
            end
        end
    end
end