function X = torow( X )
%
% TOROW Converts a vector or a matrix into a row vector.
%   If input is already a row vector, it is returned with no change.
%   If input is a column vector, it is converted into a row vector and
%   returned.
%   If input is a matrix, each column is converted into a row, and all
%   resulting rows are placed in series into a single row which is
%   returned.
%
% Input:
%   X - input vector or matrix
%
% Output:
%   X - row vector
%
% Examples:
%   torow([ 0 1 2 3 ])
%       returns [ 0 1 2 3 ]
%   torow([ 0 1 2 3 ]')
%       returns [ 0 1 2 3 ]
%   torow([ 0 1; 2 3 ])
%       returns [ 0 2 1 3 ]
%   torow([ 0 1; 2 3 ]')
%       returns [ 0 1 2 3 ]
%
% Author:	Tashi Ravach
% Version:	1.0
% Date:     07/07/2010
%

    % check if input is a vector
    [ m, n ] = size(X);
    if m == 1
        return % input is already a row vector with n columns
    elseif n==1
        X = X'; % input is converted from column vector to row vector
    elseif (m*n>n) || (m*n>m)
        X = X(:)'; % input is converted from matrix to row vector by column
    else
        X = []; % input is unknown and an empty output is returned
    end
    
end