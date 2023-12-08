function base_anim_geom = generateBaseSalp_GeneralCurvature(salp, pose, shape)
% Create base link geometry for the salp

    NThrusters = salp.geometry.thrusters.nThrusters;

    % Create base salp geometry storage
    base_anim_geom = struct();
    base_anim_geom.ticks = createBaseSalpComponent(salp, salp.visualization.colors.grid_color, 1);
    base_anim_geom.links = createBaseSalpComponent(salp, salp.visualization.colors.link_color, 1);
    base_anim_geom.thrusters = createBaseSalpComponent(salp, salp.visualization.colors.thruster_color, NThrusters);
    base_anim_geom.in_force_jets = createBaseSalpComponent(salp, salp.visualization.colors.in_jet_color, NThrusters);
    base_anim_geom.out_force_jets = createBaseSalpComponent(salp, salp.visualization.colors.out_jet_color, NThrusters);

    %% Set points for system
   
    % Base grid to draw background lines that will help judge swimmer 
    % motion
    ticks = cell(2, 1);
    plotMin = ticks;
    plotMax = plotMin;
    backgroundWidth = salp.visualization.gridSpacing;
    for idx = 1:2
        plotMin{idx} = pose(idx) - salp.geometry.length*salp.visualization.VidHalfWidth;
        plotMax{idx} = pose(idx) + salp.geometry.length*salp.visualization.VidHalfWidth;
        ticks{idx} = [-1*fliplr(0:backgroundWidth:abs(plotMin{idx})),backgroundWidth:backgroundWidth:abs(plotMax{idx})];
        
    end

    base_anim_geom.ticks.points = ticks;
    base_anim_geom.ticks.plot_min = plotMin;
    base_anim_geom.ticks.plot_max = plotMax;

    % Create overall base n-link geometry using fat chain
    jet_points = cell(NThrusters, 1);
    [B,h] = fat_backbone(salp.geometry,shape);

    base_anim_geom.links.points = B;

    thrusterCenterTransforms = cell(NThrusters,1);
    %Generate thrusters
    for idx = 1:NThrusters

        % Build thruster geometry
        thruster_radius = salp.geometry.length/NThrusters/8;
        thruster_thetas = linspace(-pi/2,pi/2,50);
        base_anim_geom.thrusters.points{idx} = [...
            thruster_radius*cos(thruster_thetas);...
            thruster_radius*sin(thruster_thetas);...
            ones(1,numel(thruster_thetas))];

        % Build force jet geometry
        thrust_jet_radius = thruster_radius*.8;
        max_thrust_length = 2/3*salp.geometry.length/NThrusters;
        thrust_jet_thetas = linspace(pi/2,3*pi/2,50);
        jet_points{idx} = [...
            max_thrust_length*cos(thrust_jet_thetas);...
            thrust_jet_radius*sin(thrust_jet_thetas);...
            ones(1,numel(thrust_jet_thetas))];
        base_anim_geom.in_force_jets.points{idx} = .5*jet_points{idx};
        base_anim_geom.in_force_jets.points{idx}(3,:) = ones(1,numel(thrust_jet_thetas));
        base_anim_geom.out_force_jets.points{idx} = jet_points{idx};

        %Store transform to thruster locations from salp body frame
        thrusterCenterTransforms{idx} = vec_to_mat_SE2(h(salp.geometry.thrusters.thrusterArcLengths(idx)));

    end
    base_anim_geom.thrusters.transforms = thrusterCenterTransforms;
    
end 