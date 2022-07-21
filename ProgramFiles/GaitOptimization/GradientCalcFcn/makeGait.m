%Makes gait definition structure of cell functions from
%fourier coefficients y for arbitrary dimension.
function gait = makeGait(y)
    ndim = size(y,2);
    
    gait.phi_def = cell(1,ndim);
    gait.dphi_def = cell(1,ndim);
    gait.ddphi_def = cell(1,ndim);
    
    for dim = 1:ndim
        gait.phi_def{dim} = @(t) phi_def(t,y(:,dim));
        gait.dphi_def{dim} = @(t) dphi_def(t,y(:,dim));
        gait.ddphi_def{dim} = @(t) ddphi_def(t,y(:,dim));
    end
                
end

function phi=phi_def(t,y)
phi = 0;
w=y(end,1);
for k = 1:1:length(y)-1
    order=floor(k/2);
    if k == 1
        phi = phi + y(k);
    elseif mod(k,2) == 0
        phi = phi + y(k)*cos(order*w*t);
    else
        phi = phi + y(k)*sin(order*w*t);
    end
end
end

function dphi=dphi_def(t,y)
dphi = 0;
w=y(end,1);
for k = 1:1:length(y)-1
    order=floor(k/2);
    if k == 1
        dphi = 0;
    elseif mod(k,2) == 0
        dphi = dphi - order*w*y(k)*sin(order*w*t);
    else
        dphi = dphi + order*w*y(k)*cos(order*w*t);
    end
end
end

function ddphi=ddphi_def(t,y)
ddphi = 0;
w=y(end,1);
for k = 1:1:length(y)-1
    order=floor(k/2);
    if k == 1
        ddphi = 0;
    elseif mod(k,2) == 0
        ddphi = ddphi - order^2*w^2*y(k)*cos(order*w*t);
    else
        ddphi = ddphi - order^2*w^2*y(k)*sin(order*w*t);
    end
end
end