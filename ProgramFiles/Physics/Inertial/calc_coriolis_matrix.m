function C = calc_coriolis_matrix(dM_alphadalpha,shape,dshape)
    % Start with (dM_alpha/dalpha*qdot)*qdot terms
    C_temp = zeros(length(shape));
    if isa(dM_alphadalpha{1},'sym')
        C_temp = sym(C_temp);
    end
    for i = 1:length(dM_alphadalpha)
        C_temp = C_temp + dM_alphadalpha{i}*dshape(i);
    end
    C = C_temp*dshape(:);
    % Add on the (-1/2)*qdot'*dM_alpha/dalpha*qdot terms
    C_temp = zeros(length(shape),1);
    if isa(dM_alphadalpha{1},'sym')
        C_temp = sym(C_temp);
    end
    for i = 1:length(dM_alphadalpha)
        C_temp(i) =  -(1/2)*dshape(:)'*dM_alphadalpha{i}*dshape(:);
    end
    C = C + C_temp;
    if isa(dM_alphadalpha{1},'sym')
        C = simplify(C);
    end
end