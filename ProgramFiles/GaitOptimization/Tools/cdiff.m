function dx = cdiff( x )
% dx = CDIFF(x) Central differencing (2nd order) of x
%   ==========
%   dx = [d/dt]x(t)
%   ==========
%    x: row of column vectors
%   ----------

[n,m] = size(x);    % find the size we must match
dx = zeros(n,m);  % preallocate

% operating on a column of row vectors
% dx(:,1)   = x(:,2)-x(:,1);     % first diff from 2 points
% dx(:,m)   = x(:,m)-x(:,m-1);   % last diff from 2 points
% dx(:,1) = nan(n,1);
% dx(:,m) = nan(n,1);

% all interior diffs from 3 points
idx = 2:m-1;
dx(:,idx) = 0.5 * (x(:,idx+1)-x(:,idx-1));

end