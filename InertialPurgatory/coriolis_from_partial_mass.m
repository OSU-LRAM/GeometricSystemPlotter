function C_alpha = coriolis_from_partial_mass(dMdalpha,jointvelocities)
    C_from_vel = zeros(size(dMdalpha));
    jointvelocities = jointvelocities(:);
    if isa(dMdalpha,'sym')
        C_from_vel = sym(C_from_vel);
    end
    C_from_pos = C_from_vel;
    for i = 1:numel(jointvelocities)
        C_from_vel = C_from_vel + dMdalpha{i}*jointvelocities(i);
        C_from_pos(i) = -(1/2)*jointvelocities'*dMdalpha{i}*jointvelocities;
    end
    C_from_vel = C_from_vel*jointvelocities;
    C_alpha = C_from_pos + C_from_vel;
    if isa(C_alpha,'sym')
        C_alpha = simplify(C_alpha,10);
    end
end