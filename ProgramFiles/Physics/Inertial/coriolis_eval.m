function C = coriolis_eval(M,dMdq,A,jointangles,jointvelocities)
    dmdalpha = pull_back_partial_mass(M,dMdq,A,jointangles);
    C_from_vel = zeros(size(dmdalpha));
    jointvelocities = jointvelocities(:);
    if isa(dmdalpha,'sym')
        C_from_vel = sym(C_from_vel);
    end
    C_from_pos = C_from_vel;
    for i = 1:numel(jointvelocities)
        C_from_vel = C_from_vel + dmdalpha{i}*jointvelocities(i);
        C_from_pos(i) = -(1/2)*jointvelocities'*dmdalpha{i}*jointvelocities;
    end
    C_from_vel = C_from_vel*jointvelocities;
    C = C_from_pos + C_from_vel;
    if isa(C,'sym')
        C = simplify(C,10);
    end
end

function dmdalpha = pull_back_partial_mass(M,dMdq,A,jointangles)
    dAdalpha = {zeros(size(A))};
    num_joints = size(A,2);
    dAdalpha = repmat(dAdalpha,1,num_joints);
    A_interp = zeros(size(A));
    for i = 1:numel(A)
        % Need some way to make this number-of-joints-agnostic
        [da1,da2] = gradient(A{i});
        dA = {da1,da2};
        for joint = 1:size(A,2)
            dAdalpha{joint}(i) = interpn(X,Y,dA{joint},jointangles(1),jointangles(2));
        end
        A_interp(i) = interpn(X,Y,A{i},jointangles(1),jointangles(2));
    end
    
    dmdalpha = cell(size(A,2),1);
    for i = 1:numel(dmdalpha)
        
        dmdalpha{i} = [-dAdalpha{i}', zeros(size(A,2))]*M*[-A_interp; eye(size(A,2))] ...
            + [-A_interp', eye(size(A,2))]*dMdq*[-A_interp; eye(size(A,2))] ...
            + [-A_interp', eye(size(A,2))]*M*[-dAdalpha{i}; zeros(size(A,2))];
    end
end