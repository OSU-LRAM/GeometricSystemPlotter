function dMdq = partial_mass_matrix(J_full,dJdq,local_inertias,base_type)
num_links = length(J_full);
num_q = size(dJdq{1},1);
% switch base_type
%     case 'fixed'
%         num_joints = num_links;
%     case 'mobile'
%         num_joints = num_links - 1;
% end

dMtemp = zeros(num_q,num_q);
% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if or(isa(J_full{1},'sym'),isa(local_inertias{1},'sym'))
    dMtemp = sym(dMtemp);
    cell2func = @cell2sym;
else
    cell2func = @cell2mat;
end
dMdq = cell(num_q,1); dMdq(1:end) = {dMtemp};
for q = 1:num_q
    for link = 1:num_links
        dJtemp = cell2func(dJdq{link}(:,q)');
         dMdq{q} = dMdq{q} + dJtemp' * local_inertias{link} * J_full{link} ...
            + J_full{link}' * local_inertias{link} * dJtemp;
    end
    if isa(dMtemp,'sym')
        dMdq{q} = simplify(dMdq{q},'Steps',10);
    end
end