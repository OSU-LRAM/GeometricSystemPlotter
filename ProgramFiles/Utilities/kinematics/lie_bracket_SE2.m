function lie_bracket = lie_bracket_SE2(a,b)
% Calculates the lie bracket of two vectors a, b which are members of the
% set SE2. The vectors a and b should be 3x1 vectors of the form:
%       a = [a_x; a_y; a_theta]
%       b = [b_x; b_y; b_theta]
% with the subscripts denoting the x, y, and theta coordinates of the
% vectors, respectively.

% The lie bracket of two vectors finds the net effect of making sequential,
% infinitesimal moves in the a, b, -a, -b directions; it returns a vector 
% whose values are obtained by using a small angle approximation of TeLg 
% around theta = 0. This is equivalent to the conventional lie bracket
% [X,Y] = XY - YX where X and Y are SE2  matrices with small angle
% approximations used.

lie_bracket = [b(3)*a(2) - a(3)*b(2);   % x-direction term
               a(3)*b(1) - b(3)*a(1);   % y-direction term
               0];                      % theta-direction term
end