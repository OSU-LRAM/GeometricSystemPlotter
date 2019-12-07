% TODO: find out how to un-hardcode this
%home = "/home/qkonyn/Programming_Projects/GeometricSystemPlotter/";
%addpath(home + "UserFiles/GenericUser/Systems");
%addpath(home + "ProgramFiles/Geometry");

s = sysf_two_link_lowRe;

friction_direction = 1; % todo: vary this
shape_step = 0.5;
% shape is alpha; shapechange is alpha dot
for shape = s.grid_range(1):shape_step:s.grid_range(2)
    for shapechange = -0.5:0.1:0.5
        % get A(alpha)
        [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry, s.physics, shape);

        % gcirc right
        body_vel = A * shapechange;

        % get Jfull(alpha)
        [~, ~, J_full, ~, ~] = N_link_chain(s.geometry, shape);

        for link = 1:length(s.geometry.linklengths)
            body_vel = J_full{link} * [body_vel; shapechange];
            dx_direction = sign(body_vel(1)); %index 1 is x, right?
            agreement = friction_direction * dx_direction;
            % todo: store result for plotting
        end
    end
end