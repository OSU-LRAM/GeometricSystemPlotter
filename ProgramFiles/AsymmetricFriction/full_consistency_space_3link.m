function [consistent_system,consistency_count,consistent_A] = full_consistency_space_3link(s, adot_scale)
%Returns a 4D space of which Smooth/Rough subsystem the scaled swimmer is
%consistent with at each combination of alphas and alpha dots
% directions_key = ["SSS","SSR","SRS","SRR","RSS","RSR","RRS","RRR","AXES"];

nlinks = length(s.geometry.linklengths);

a_grid = s.grid.eval{1};
% adot_scale = 0.1;
adot_grid = s.grid.eval{1}*adot_scale;

directions_possible = ff2n(nlinks);

consistent_system = zeros(length(a_grid), length(a_grid), length(adot_grid), length(adot_grid));
consistency_count = zeros(length(a_grid), length(a_grid), length(adot_grid), length(adot_grid));
consistent_A = zeros(length(a_grid), length(a_grid), length(adot_grid), length(adot_grid), 3, 2);

% loop through shapes
for j = 1:size(s.grid.eval{1},1)
    for k = 1:size(s.grid.eval{1},2)
        shape = [s.grid.eval{1}(j,k);s.grid.eval{2}(j,k)];
        [~, ~, J_full, ~, ~] = N_link_chain(s.geometry, shape);
        % loop through every combination of Rough/Smooth links
        for i = 1:size(directions_possible,1)
            links_backwards = directions_possible(i,:);
            friction_directions = ones(1,nlinks) - 2 * links_backwards;
            [A, ~, ~, ~, ~] = LowRE_connection_discrete(s.geometry,s.physics, shape, links_backwards);
            % loop through shapechanges
            for l = 1:size(s.grid.eval{1},1)
                for m = 1:size(s.grid.eval{1},2)
                    shapechange = adot_scale*[s.grid.eval{1}(l,m);s.grid.eval{2}(l,m)];
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
                        consistent_system(j, k, l, m) = i;
                        consistent_A(j, k, l, m, :, :) = A;
                        consistency_count(j, k, l, m) = consistency_count(j, k, l, m) + 1;
                    end
                    if consistency_count(j, k, l, m) == 8
                        % special code 9 for when all systems are
                        % consistent (axes)
                        consistent_system(j, k, l, m) = 9;
                    end
                end
            end
    end
    end
end

end

