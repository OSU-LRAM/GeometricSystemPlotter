function y = path_from_fourier(f,n,dimension)
% Returns the shape space parametrization of the gait at n points when provided with
% the fourier coefficients f. The gait is returned as a self-closed gait
% (i.e. the first and last rows of the output y are the same point).
% Inputs:
%   f: Fourier coefficients that parametrize the gait.
%   n: Number of points that should compose the gait, less the final
%       self-closed point
%   dimension: Number of shape variables for system

    y = zeros(n+1,dimension);
    % Determine time period based on value of fourier frequency
    w = f(end,1);
    T = 2*pi/w;
    % Create time vector at which to evaluate points of gait
    t = linspace(0,T,n+1);
    % Evaluate the shape-space parametrization of the gait at every time
    % value in t
    for j=1:dimension
        for i=1:1:n+1
            y(i,j)=f(1,j)+f(2,j)*cos(w*t(i))+f(3,j)*sin(w*t(i))+f(4,j)*cos(2*w*t(i))+...
                +f(5,j)*sin(2*w*t(i))+f(6,j)*cos(3*w*t(i))+f(7,j)*sin(3*w*t(i))+...
                +f(8,j)*cos(4*w*t(i))+f(9,j)*sin(4*w*t(i));
        end
    end
end