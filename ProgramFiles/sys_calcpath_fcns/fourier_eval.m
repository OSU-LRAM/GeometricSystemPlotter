function f = fourier_eval(x,AN,BN,T)
% Evaluate a fourier function with coefficients AN and BN at points x for a
% period T

ii = (1:length(AN))-1;   % The indices. (with zero-indexing)
lgh = length(x);
[ii x] = meshgrid(ii,x); % Arrays for angles.
AN = repmat(AN',lgh,1);   % Create arrays for approximation.
BN = repmat(BN',lgh,1);
theta = ii.*x.*(2*pi)/T; % The angles.
f = (sum(AN.*cos(theta)+BN.*sin(theta),2))'; % the F approximation.

end