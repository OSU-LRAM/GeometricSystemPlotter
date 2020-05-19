function dJdq = mobile_jacobian_derivative(J_full)
% Calculates the partial of the full Jacobian J_full with respect to the
% configuration variables q = [g_circ r]', where g_circ is the system's
% global coordinate configuration and r is the local shape-space
% configuration.
%
% Inputs for a system with m links and q of size (n by 1):
%
%   J_full: Full Jacobian matrix for the system, including entries for all
%       configuration variables q; (m by 1) cell array where each cell contains
%       a (3 by n) matrix representing the Jacobian of link m.
%
% Output:
%
%   dJdq: Partial of full Jacobian with respect to the configuration
%       variables q; (m by 1) cell array where each entry corresponds to
%       the partial derivative of the Jacobian for link m, which is another
%       cell array of size (n by n) where each entry (j,k) corresponds to a
%       three-element vector that represents the partial of Jacobian column
%       j with respect to configuration variable k.

% Determine how many joint angles are in the system; the first three terms
% in the Jacobian are for the body coordinates, with any additional columns
% corresponding to the joint space of the system
n_joints = size(J_full{1},2) - 3;
m_links = length(J_full);

% Prototype a template to fit into each level of the output cell array
dJdq = cell(m_links,1);
template = cell(3+n_joints,3+n_joints);
zero_vec = zeros(3,1);
% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if isa(J_full{1},'sym')
    zero_vec = sym(zero_vec);
end
template(1:end) = {zero_vec};

dJdq(1:end) = {template};


for link = 1:m_links
% The Jacobian derivative dJ/dq is of the form
%
%       dJ/dq = [dJg/dRg, dJg/dalpha; dJalpha/dRg, dJalpha/dalpha]
%
% where dJg/dRg is a 3x3 augmented matrix of three-element zero vectors,
% dJg/dalpha is a 3xm augmented matrix with each element i,j as the lie 
% bracket [Jg_i,Jalpha_j] for each body coordinate i and joint j,
% dJalpha/dRg is a mx3 augmented matrix of three-element zero vectors, and
% dJalpha/dalpha is the lie bracket [Jalpha_i,Jalpha_j] for each joint i
% and j with j > i.

    % No calculation necessary for dJg/dRg and dJalpha/dRg - leave dJdq in
    % these areas as zero vectors

    % Calculating dJg/dalpha
    for idx = 1:3
        for idx2 = 1:n_joints
            dJdq{link}{idx,3+idx2} = lie_bracket_SE2(J_full{link}(:,idx),J_full{link}(:,3+idx2));
        end
    end

    % Calculating dJalpha/dalpha
    for idx = 1:n_joints
        for idx2 = idx+1:n_joints
            dJdq{link}{3+idx,3+idx2} = lie_bracket_SE2(J_full{link}(:,3+idx),J_full{link}(:,3+idx2));
        end
    end
end
end