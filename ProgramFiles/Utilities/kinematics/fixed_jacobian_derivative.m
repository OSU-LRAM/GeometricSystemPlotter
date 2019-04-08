function dJdalpha = fixed_jacobian_derivative(J,jointvelocities)
% Determine how many links there are in the Jacobian
m_links = length(J);
% Find how many joints are in the system
n_joints = size(J{1},2);
% Makes sure jointvelocities is a column vector
jointvelocities = jointvelocities(:);

% The Jacbian will have a derivative term for each link and joint.
% Preallocate an array to hold the resulting derivative, which will be the
% same size as the Jacobian
dJdalpha = cell(m_links,1);
dJdalpha(1:m_links) = {zeros(size(J{1}))};
% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if or(isa(J,'sym'),isa(jointvelocities,'sym'))
    dJdalpha = sym(dJdalpha);
end

% Calculate the Jacobian for each link
for link = 1:length(J)
    % Preallocate a temporary array to hold the resulting lie brackets
    lie_temp = zeros(3*n_joints,n_joints);
    % If we're working with symbolic variables, then we need to explicitly make
    % the array symbolic, because matlab tries to cast items being inserted
    % into an array into the array class, rather than converting the array to
    % accomodate the class of the items being inserted 
    if or(isa(J,'sym'),isa(jointvelocities,'sym'))
        lie_temp = sym(lie_temp);
    end

    % For the derivative, only links down the chain from the joint
    % currently being looked at (i.e., the current row) are effected, so
    % only calculate the lie bracket for the upper triangle
    for row = 1:1:n_joints
        for col = row+1:1:n_joints
            lie_temp(3*row-2:3*row,col) = lie_bracket_SE2(J{link}(:,row),J{link}(:,col));
        end
    end
    
    % Multiply the lie bracket result vectors by the joint velocities and
    % return the Jacobian derivative by reshaping the resulting vector to
    % be of the same shape as the Jacobian.
    dJtemp = lie_temp*jointvelocities;
    dJdalpha{link} = reshape(dJtemp,3,n_joints);
end

end % jacobian_derivative