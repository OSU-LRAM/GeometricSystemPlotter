function dMdq = partial_mass_matrix(J,dJdq,local_inertias)
q_num = length(J);
dMtemp = zeros(q_num,q_num);
% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if or(isa(J{1},'sym'),isa(local_inertias{1},'sym'))
    dMtemp = sym(dMtemp);
    cell2func = @cell2sym;
else
    cell2func = @cell2mat;
end
dMdq = cell(q_num); dMdq(1:end) = {dMtemp};
for q = 1:q_num
    for link = 1:q_num
        dJtemp = cell2func(dJdq{link}(:,q)');
        dMdq{q} = dMdq{q} + dJtemp' * local_inertias{link} * J{link} ...
            + J{link}' * local_inertias{link} * dJtemp;
    end
    if isa(dMtemp,'sym')
        dMdq{q} = simplify(dMdq{q},'Steps',10);
    end
end
end