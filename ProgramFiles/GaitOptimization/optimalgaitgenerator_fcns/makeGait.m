%Makes gait definition structure of cell functions from
%fourier coefficients y for arbitrary dimension.
function gait = makeGait(y)

    ndim = size(y,2);
    
    gait.phi_def = cell(1,ndim);
    gait.dphi_def = cell(1,ndim);
    gait.ddphi_def = cell(1,ndim);
    
    for dim = 1:ndim
        
        w = y(end,dim);
        gait.phi_def{dim} = @(t) y(1,dim)+y(2,dim)*cos(w*t)+y(3,dim)*sin(w*t)+y(4,dim)*cos(2*w*t)+...
                                +y(5,dim)*sin(2*w*t)+y(6,dim)*cos(3*w*t)+y(7,dim)*sin(3*w*t)+...
                                +y(8,dim)*cos(4*w*t)+y(9,dim)*sin(4*w*t);
        gait.dphi_def{dim} = @(t) -w*y(2,dim)*sin(w*t)+w*y(3,dim)*cos(w*t)-2*w*y(4,dim)*sin(2*w*t)+...
                                  +2*w*y(5,dim)*cos(2*w*t)-3*w*y(6,dim)*sin(3*w*t)+3*w*y(7,dim)*cos(3*w*t)+...
                                  -4*w*y(8,dim)*sin(4*w*t)+4*w*y(9,dim)*cos(4*w*t);
        gait.ddphi_def{dim} = @(t) -w^2*y(2,dim)*cos(w*t)-w^2*y(3,dim)*sin(w*t)-4*w^2*y(4,dim)*cos(2*w*t)+...
                                   -4*w^2*y(5,dim)*sin(2*w*t)-9*w^2*y(6,dim)*cos(3*w*t)-9*w^2*y(7,dim)*sin(3*w*t)+...
                                   -16*w^2*y(8,dim)*cos(4*w*t)-16*w^2*y(9,dim)*sin(4*w*t);
        
    end
                
end