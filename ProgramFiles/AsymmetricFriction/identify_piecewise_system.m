function [consistent_system, consistency_count] = identify_piecewise_system(s, a, adot)
%IDENTIFY_PIECEWISE_SYSTEM Summary of this function goes here
%   Find the system that's consistent at each point in the alpha / alphadot
%   space.
%   TODO: If there's any situation where no system is consistent, or
%   multiple are consistent, throw an error.
%   Output:
%       0: no system is consistent
%       1: first system is consistent
%       ...
%       last system: any system is consistent... occurs on axes

consistent_system = zeros(length(a), length(adot));
consistency_count = zeros(length(a), length(adot)); % for error checking

for dir = 0:3
    backwards = [mod(dir,2), floor(dir/2)]; % this is a hacky way to do all 4 direction combos with binary
    friction_direction = [1, 1] - 2 * backwards;
    for i = 1:length(a)
        shape = a(i);
        
        % todo: the also don't depend on dir, so make this the outer loop
        % these two don't depend on alpha dot:
        % get A(alpha)
        [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry, s.physics, shape, backwards);
        % get Jfull(alpha)
        [~, ~, J_full, ~, ~] = N_link_chain(s.geometry, shape);
        
        for j = 1:length(adot)
            shapechange = adot(j);
            
            % gcirc right
            body_vel = A * shapechange;
            
            all_links_agree = 1;
            for link = 1:length(s.geometry.linklengths)
                link_body_vel = J_full{link} * [body_vel; shapechange];
                dx_direction = sign(link_body_vel(1));
                if(dx_direction == 0) % non-if way to do this?
                    agreement = 1;
                else
                    agreement = friction_direction(link) * dx_direction;
                end
                all_links_agree = all_links_agree && (1 + agreement);
            end
            if all_links_agree
                consistent_system(i, j) = 1 + dir;
                consistency_count(i, j) = consistency_count(i,j) + 1; % no increment operator?!?!?!?
            end
        end
    end
end

end

