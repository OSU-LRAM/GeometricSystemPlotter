function A = get_correct_A(s, shape, shapechange)
%GET_CORRECT_A Summary of this function goes here
%   Detailed explanation goes here

nlinks = length(s.geometry.linklengths);

if nlinks == 1
    % RANGES BETWEEN -PI AND PI ONLY
    opening = sign(shape) == sign(shapechange);

    if opening
        dir = 3 - 1;
    else
        dir = 2 - 1;
    end

    backwards = [mod(dir,2), floor(dir/2)]; % this is a hacky way to do all 4 direction combos with binary

    % get A(alpha)
    [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry, s.physics, shape, backwards);
else % 2 links
    if s.physics.drag_bw_ratio == 1
        [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry,s.physics, shape);
    else
        directions_possible = ff2n(nlinks);
        [~, ~, J_full, ~, ~] = N_link_chain(s.geometry, shape);
        % go through each possible combo of directions to find a consistent one
        for i = 1:size(directions_possible,1)
            links_backwards = directions_possible(i,:);
            friction_directions = ones(1,nlinks) - 2 * links_backwards;
            [potential_A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry,s.physics, shape, links_backwards);
            body_vel = -potential_A * shapechange;
            all_links_agree = 1;
            for link = 1:length(s.geometry.linklengths)
                link_body_vel = J_full{link} * [body_vel; shapechange];
                dx_direction = sign(link_body_vel(1));
                if(dx_direction == 0) % non-if way to do this?
                    agreement = 1;
                else
                    agreement = friction_directions(link) * dx_direction;
                end
                all_links_agree = all_links_agree && (1 + agreement);
            end
            if all_links_agree
                A = potential_A;
            end
        end
    end
    if ~exist('A')
        shape
        shapechange
        J_full
    end
end

