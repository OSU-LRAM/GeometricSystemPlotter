function C = calc_coriolis_matrix(dM,dshape)
    % Start with (dM_alpha/dalpha*qdot)*qdot terms
    C_temp = zeros(length(dshape));
    if isa(dM{1},'sym')
        C_temp = sym(C_temp);
    end
    for i = 1:length(dM)
        C_temp = C_temp + dM{i}*dshape(i);
    end
    C = C_temp*dshape(:);
    % Add on the (-1/2)*qdot'*dM_alpha/dalpha*qdot terms
    C_temp = zeros(length(dshape),1);
    if isa(dM{1},'sym')
        C_temp = sym(C_temp);
    end
    for i = 1:length(dM)
        C_temp(i) =  -(1/2)*dshape(:)'*dM{i}*dshape(:);
    end
    C = C + C_temp;
    if isa(dM{1},'sym')
        C = simplify(C);
    end
end