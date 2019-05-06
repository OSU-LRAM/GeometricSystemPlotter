function dJdq = fixed_jacobian_derivative(J)
% Determine how many joint angles are in the system; the first three terms
% in the Jacobian are for the body coordinates, with any additional columns
% corresponding to the joint space of the system
n_joints = size(J{1},2);
m_links = length(J);

dJdq = cell(m_links,1);
template = cell(n_joints,n_joints);
zero_vec = zeros(3,1);
% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if isa(J{1},'sym')
    zero_vec = sym(zero_vec);
end
template(1:end) = {zero_vec};
dJdq(1:end) = {template};


% Calculate the Jacobian for each link
for link = 1:m_links
    % Calculating dJalpha/dalpha
    for idx = 1:n_joints
        for idx2 = idx+1:n_joints
            dJdq{link}{idx,idx2} = lie_bracket_SE2(J{link}(:,idx),J{link}(:,idx2));
        end
    end
end

end % jacobian_derivative