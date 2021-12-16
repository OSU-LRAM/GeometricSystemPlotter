function [consistent_system, consistent_A] = determine_subsystem(s,shape,shapechange,J_full)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nlinks = length(s.geometry.linklengths);
directions_possible = ff2n(nlinks);

consistency_count = 0;

for i = 1:size(directions_possible,1)
    links_backwards = directions_possible(i,:);
    friction_directions = ones(1,nlinks) - 2 * links_backwards;
    [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry,s.physics, shape, links_backwards);
    
    body_vel = -A * shapechange;

    all_links_agree = 1;
    for link = 1:length(s.geometry.linklengths)
        link_body_vel = J_full{link} * [body_vel; shapechange];
        dx_direction = sign(link_body_vel(1));
        if(dx_direction == 0) % non-if way to do this?
            agreement = 1; % if the link isn't moving, just say it agrees either way
        else
            agreement = friction_directions(link) * dx_direction;
        end
        all_links_agree = all_links_agree && (1 + agreement);
    end
    if all_links_agree
        consistent_system = i;
        consistent_A = A;
        consistency_count = consistency_count + 1;
    end
    if consistency_count == 8
        % special code 9 for when all systems are
        % consistent (axes)
        consistent_system = 9;
    end
end

end

