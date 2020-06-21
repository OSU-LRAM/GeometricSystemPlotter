function ddJdqq = mobile_jacobian_second_derivative(J_full,dJdq)
% Determine how many joint angles are in the system; the first three terms
% in the Jacobian are for the body coordinates, with any additional columns
% corresponding to the joint space of the system
n_joints = size(J_full{1},2) - 3;
m_links = length(J_full);

% We need ddJdqq for every link, and each link's second derivative will be
% in terms of dRg and dalpha.
ddJdqq = cell(m_links,1);
template = cell(3+n_joints,3+n_joints,3+n_joints);
template(1:end) = {zeros(3,1)};
ddJdqq(1:end) = {template};
% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if isa(J_full{1},'sym')
    ddJdqq = sym(ddJdqq);
end

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
% The second derivative of the Jacobian is simply dJ/dq with additional
% partial derivatives d/dRg and d/dalpha applied. For notation purposes,
% these partials are kept as separate sheets in the third dimension of
% ddJdqq.

    % No calculation is necessary for the second partial for the d/dRg
    % sheets, as the partials of these terms is always zero. Leave these 
    % (the first three sheets) as the zero vectors from the template.

    % The partial terms dJg/dRg and dJalpha/dRg are always zero; thus no
    % calculation is necessary for the second partial d/dalpha in the first
    % three columns of the fourth sheet and onward. Leave these as zero 
    % vectors from the template.

    % Calculating ddJg/dalpha^2
    for idx = 1:3 % Each of Jg1, Jg2, Jg3
        for idx2 = 1:n_joints % Each dJg/dalpha
            for idx3 = 1:n_joints % Second d/dalpha for each alpha
                ddJdqq{link}{idx,3+idx2,3+idx3} = ...
                    lie_bracket_SE2(dJdq{link}{idx,3+idx2},J_full{link}(:,3+idx3));
            end
        end
    end

    % Calculating ddJalpha/dalpha^2
    for idx = 1:n_joints % Each of Jalpha1, Jalpha2, ...
        for idx2 = idx+1:n_joints % Each dJalpha/dalpha 
            for idx3 = idx+1:n_joints % Second d/dalpha for each alpha
                ddJdqq{link}{3+idx,3+idx2,3+idx3} = ...
                    lie_bracket_SE2(dJdq{link}{3+idx,3+idx2},J_full{link}(:,3+idx3));
            end
        end
    end
end
end