function ddx = ddiff( x )
% x: row of column vectors

[n,m] = size(x);    % find the size we must match
ddx = nan(n,m);  % preallocate

ddx(:,2:end-1) = x(:,3:end) - 2*x(:,2:end-1) + x(:,1:end-2);

end

