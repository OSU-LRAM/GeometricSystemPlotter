% find_g_sys.m
% sysplotter wrapper for integrated position, respecting config. space
% note: SO(3) positions are in matrix form, pending choice of coordinate

function p = find_g_sys(s,p)
    % switch between integrators based on group
    if s.conf_space == LieGroups.SE2
        % use working legacy setup for now
        % this ought to be consolidated during rewrite
        p = find_g(s,p);
    elseif s.conf_space == LieGroups.SO3
        % redundant so3 wrapper (below) for simplicity
        p = find_g_so3(s,p);
    end
end

% so3 integration wrapper (produces same outputs as find_g)
function p = find_g_so3(s,p)
    % respect multiple paths
    for i = 1:length(p.phi_fun)
        % integrate segments separately
        for j = 1:length(p.phi_fun{i})
            % get local g with so3 fixed-step integrator
            [~, g_all] = so3_integrator(p.time{i}{j},...
                                        p.phi_fun{i}{j},...
                                        p.dphi_fun{i}{j},...
                                        s.grid.eval,...
                                        s.vecfield.eval.content.Avec);
            [~, g_all_opt] = so3_integrator(p.time{i}{j},...
                                            p.phi_fun{i}{j},...
                                            p.dphi_fun{i}{j},...
                                            s.grid.eval,...
                                            s.vecfield.eval.content.Avec_optimized);
            p.G_locus{i,j}.G_local = g_all;
            p.G_locus{i,j}.G_opt_local = g_all_opt;
            % lift each point into previous frame
            if j == 1
                % no previous frame
                p.G_locus{i,j}.G = p.G_locus{i,j}.G_local; %init. cond.?
                p.G_locus{i,j}.G_opt = p.G_locus{i,j}.G_opt_local;
                continue;
            end
            % multiply ea. point by end of last segment
            % this notation is awful, can we improve with another pass?
            p.G_locus{i,j}.G = cell(size(p.G_locus{i,j}.G_local));
            p.G_locus{i,j}.G_opt = cell(size(p.G_locus{i,j}.G_opt_local));
            for k = 1:length(p.G_locus{i,j-1})
                p.G_locus{i,j}.G{k} = p.G_locus{i,j-1}.G{end} * p.G_locus{i,j}.G_local{k};
                p.G_locus{i,j}.G_opt{k} = p.G_locus{i,j-1}.G_opt{end} * p.G_locus{i,j}.G_opt_local{k};
            end
        end
        % concat all configurations
        p.G_locus_full{i}.G = [];
        p.G_locus_full{i}.G_opt = [];
        for j = 1:length(p.phi_fun{i})
            p.G_locus_full{i}.G = [p.G_locus_full{i}.G, p.G_locus{i,j}.G];
            p.G_locus_full{i}.G_opt = [p.G_locus_full{i}.G_opt, p.G_locus{i,j}.G_opt];
        end
        % compute cBVI
        % should be general-purpose already
        if ~strcmp(p.cBVI_method{i}{1},'none')
            [p.cBVI, p.cBVI_opt] = integrate_cBVI(s,p);
        end
    end
end