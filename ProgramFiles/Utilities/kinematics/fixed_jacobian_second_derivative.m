function ddJdqq = fixed_jacobian_second_derivative(J,dJdq)

n_joints = size(J{1},2);
m_links = length(J);

ddJdqq = cell(m_links,1);
template = cell(n_joints,n_joints,n_joints);
zero_vec = zeros(3,1);
% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if isa(J{1},'sym')
    zero_vec = sym(zero_vec);
end
template(1:end) = {zero_vec};
ddJdqq(1:end) = {template};

%
for link = 1:m_links
    for idx = 1:n_joints
        for idx2 = idx+1:n_joints
            for idx3 = idx+1:n_joints
                ddJdqq{link}{idx,idx2,idx3} = lie_bracket_SE2(dJdq{link}{idx,idx2},J{link}(:,idx3));
            end
        end
    end
end