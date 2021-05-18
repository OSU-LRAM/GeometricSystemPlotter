function varargout = ndgradient(varargin)
% Wrapper on gradient to work with ndgrid instead of meshgrid

% For now we assume that there is a grid being provided. If making this
% handle all of gradient's edge cases, will need to generalize
varargout = cell(1,numel(varargin)-1);

grad_permute = [2,1,3:numel(varargout)];

[varargout{grad_permute}] = gradient(varargin{1},varargin{grad_permute+1});

end