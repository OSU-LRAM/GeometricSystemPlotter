function animate_locomotor_race(export,info_needed)
    
    % Look up the geometry specification for the systems:
    for idx_system = 1:numel(info_needed.current_system2)
        sysfile = fullfile(info_needed.datapath, ['sysf_', info_needed.current_system2{idx_system}, '_calc.mat']);
        load(sysfile,'s')
        info_needed.s(idx_system) = s;
    end
        
  	% Declare a directory name file names for the movies
	destination_root = fullfile(info_needed.datapath,'Animation','swimrace');
	destination_list = fieldnames(info_needed.Movies);
            
        if strcmp(info_needed.Coordinates,'minperturbation_coords')
            destinationsuffix = 'min_perturb';
        else
            destinationsuffix = 'original_coords';
        end

	destination = cellfun(@(x) fullfile(destination_root,['race' '__' 'race' '__' x '__' destinationsuffix]),destination_list,'UniformOutput',false); 

    % Flag movies that should be generated, don't skip any movies
    export_list = cellfun(@(x) info_needed.Movies.(x), destination_list);
    skip_list = zeros(size(export_list));
    info_needed.export_list = export_list;
    
    %%%%%%%%
	% Create animation elements, and store them in the frame_info structure
	frame_info{1} = create_elements(info_needed); %setup function, defined below in file
	
	% Designate animation function
	frame_gen_function...
		= @execute_gait; % frame function, defined below in file

	% Declare timing
	timing.duration = info_needed.Duration; % in seconds
	timing.fps = info_needed.Framerate;     % create frames for specified frame rate
	timing.pacing = @(y) softspace(0,1,y);  % Use a soft start and end, using the included softspace function

    
    % Pass the system and gait names to the frame_info struct
    frame_info{1}.sysname = info_needed.current_system2;
    frame_info{1}.pathname = info_needed.current_shch2;
    frame_info{1}.s = info_needed.s;
    frame_info{1}.drawing_baseframe_inserted = repmat({0},[numel(info_needed.current_system2),1]);


	% Animate the movie
	[frame_info, endframe]...
		= sysplotter_animation_race(frame_gen_function,frame_info,timing,destination,export_list,skip_list,0);

end

function h = create_elements(info_needed)

    

	h.f = figure(17);                            % Designate a figure for this animation
	clf(h.f)                                     % Clear this figure
	
    set(h.f,...
        'position',...
        [600 1500 800 350],...
        'paperposition',[0 0 10 4.5],...
        'papersize',[10 4.5],...
        'Color','w',...
        'renderer','painters');

%     set(h.f,...
%         'position',[600 1500 800 600],...
%         'paperposition',[0 0 8 6],...
%         'papersize',[8 6],...
%         'Color','w',...
%         'renderer','painters');

	h.ax = axes('Parent',h.f);                      % Create axes for the plot
 	axis(h.ax,'equal','off');                       % Make the plot isometric and make the axes invisible
  	set(h.ax,...
        'XLim',[-1,4]*.7/.45*info_needed.s(1).geometry.length,...  % Axes scaled to system scale
        'YLim',[-2.5,2.5]*.7*info_needed.s(1).geometry.length);
 	set(h.ax,'Position',[0 0 1 1])                  % Make the axes fill the whole window

    data_source = info_needed.datapath;
	system_name = info_needed.current_system2;
	gait_name = info_needed.current_shch2;

    for idx_system = 1:numel(info_needed.current_system2)
        
        %Create the locomotor
        info_needed2 = info_needed;
        info_needed2.s = info_needed.s(idx_system);
        h.robot{idx_system} = create_locomotor(...
            h.ax,...
            data_source,...
            system_name{idx_system},...
            info_needed2);
	
        %%%%%%%%
        % Create the tracer objects
        load('sysplotter_config','Colorset')
        h.tld{idx_system} = line('Parent',h.ax,'XData',[],'YData',[],'ZData',[],'LineWidth',5,'Color',Colorset.spot,'Marker','o','MarkerFaceColor',Colorset.spot,'MarkerSize',5,'LineStyle','none');
        h.tl{idx_system} = line('Parent',h.ax,'XData',[],'YData',[],'ZData',[],'LineWidth',5,'Color',Colorset.spot); %Tracer line for translation

        % Extract the position, shape, and cost data
        load(fullfile(data_source,['sysf_' system_name{idx_system} '__shchf_' gait_name{idx_system}]));
        shaperaw{idx_system} = p.phi_locus_full{1}.shape;

        if strcmp(info_needed.Coordinates,'minperturbation_coords')
            posraw{idx_system} = p.G_locus_full{1}.G_opt;
            h.drawing_baseframe{idx_system} = info_needed.current_system2{idx_system};
        else
            posraw{idx_system} = p.G_locus_full{1}.G;
        end
        
        costraw{idx_system} = p.G_locus_full{1}.S;
        
        speedraw{idx_system} = abs(posraw{idx_system}(end,1))/costraw{idx_system}(end);
        
    end
    
    [maxspeed,fastestsystem] = max(cell2mat(speedraw));
    
    for idx_system = 1:numel(info_needed.current_system2)
        
       costscaled{idx_system} = costraw{idx_system} * maxspeed;
       speedscaled{idx_system} = speedraw{idx_system}/maxspeed;
        
    end

	
	
	% Specify the number of times to repeat the cycle, and the number of
	% negative cycle-displacements to start from
    maxdisp = 3;
	n_gaits_max = ceil(maxdisp/abs(posraw{fastestsystem}(end,1)));
    t_max = costscaled{fastestsystem}(end)*n_gaits_max;
	start_pos = 0;
	
	% get full configuration history
	[h.shapedata, h.posdata] = gait_concatenator(shaperaw,posraw,n_gaits_max,start_pos,t_max,costscaled);
	

end



%%%%%%%%%%%%%%%%%%
% Frame content to draw during the gait
function frame_info = execute_gait(frame_info,tau)

	% Call out the frame_info pieces that contain robot verse shapespace
	% data
	irobot = 1;

	
	%%%%%%
    % Get the system position and shape at fractional time tau
    for idx_system = 1:numel(frame_info{irobot}.shapedata)
    
        % Generate a core timing vector from zero to 1, with as many entries as
        % data points in the kinematic history
        timing_base = linspace(0,1,size(frame_info{irobot}.shapedata{idx_system},1));

        % Interpolate the shape variables at the specified time
        config.shape = zeros(size(frame_info{irobot}.shapedata{idx_system},2),1);
        for idx = 1:numel(config.shape)
            config.shape(idx) = interp1(timing_base,frame_info{irobot}.shapedata{idx_system}(:,idx),tau);
        end
    
        % Interpolate the position variables at the specified time
        config.x = interp1(timing_base,frame_info{irobot}.posdata{idx_system}(:,1),tau);
        config.y = interp1(timing_base,frame_info{irobot}.posdata{idx_system}(:,2),tau);
        config.theta = interp1(timing_base,frame_info{irobot}.posdata{idx_system}(:,3),tau);
        
        % offset the y values to put the systems next to each other
        yoffset = + (idx_system - 1) - (numel(frame_info{irobot}.shapedata)-1)/2;
        config.y = config.y + yoffset;
	
        %%%%
        % Place the locomotor at the position

        % If a drawing baseframe has been specified, append the original
        % system baseframe with the one that has been specified
        if isfield(frame_info{irobot},'drawing_baseframe') && ~frame_info{irobot}.drawing_baseframe_inserted{idx_system}
            if ~iscell(frame_info{irobot}.s(idx_system).geometry.baseframe)
                frame_info{irobot}.s(idx_system).geometry.baseframe ...
                    = {frame_info{irobot}.s(idx_system).geometry.baseframe};
            end
            if ~iscell(frame_info{irobot}.s(idx_system).geometry.baseframe)
                frame_info{irobot}.drawing_baseframe ...
                    = {frame_info{irobot}.drawing_baseframe};
            end
            frame_info{irobot}.s(idx_system).geometry.baseframe = [frame_info{irobot}.s(idx_system).geometry.baseframe frame_info{irobot}.drawing_baseframe{idx_system}];
            frame_info{irobot}.drawing_baseframe_inserted{idx_system} = 1;
        end
    
        % Use the configuration to place the locomotor
        frame_info{irobot}.robot{idx_system} = place_locomotor(frame_info{irobot}.robot{idx_system},config,frame_info{irobot}.s(idx_system));

        % draw the locomotor
        frame_info{irobot}.robot{idx_system} = draw_locomotor(frame_info{irobot}.robot{idx_system},0);

        % Draw the tracer dot
        set(frame_info{irobot}.tld{idx_system},'XData',[frame_info{irobot}.posdata{idx_system}(1,1) config.x],'YData',[frame_info{irobot}.posdata{idx_system}(1,2)+yoffset config.y],'ZData',[5,5]);

        % Draw the tracer
        set(frame_info{irobot}.tl{idx_system},'XData',frame_info{irobot}.posdata{idx_system}(1:floor(tau*end),1),'YData',frame_info{irobot}.posdata{idx_system}(1:floor(tau*end),2)+yoffset,'ZData',5*ones(size(frame_info{irobot}.posdata{idx_system}(1:floor(tau*end),1))));


    end
    
	% Declare a print method (in this case, print 150dpi png files of
	frame_info{irobot}.printmethod = @(dest) print(frame_info{irobot}.f,'-r 150','-painters',dest);
	

end

%%%%%%
function [shapedata, posdata] = gait_concatenator(shaperaw,posraw,n_gaits_max,start_pos,t_max,costscaled)


	
	% Store the shape and data to the frameinfo structure
    for idx_system = 1:numel(shaperaw)
        
        n_gaits = ceil(t_max/costscaled{idx_system}(end));
        
        shapedata{idx_system} = cat(1,shaperaw{idx_system}(1,:),repmat(shaperaw{idx_system}(2:end,:),[n_gaits,1]));
    
        
        %%%%%%
        % Concatenate the gait displacements together.

        % Get the displacement over one cycle
        cyclic_displacement = posraw{idx_system}(end,:);
	
        % convert this cyclic displacement into an SE(2) matrix
        cyclic_displacement_m = vec_to_mat_SE2(cyclic_displacement);

        % Add a set of ones onto the xy values, so that they can be
        % translated/rotated by SE(2) matrices
        posraw_xy_augmented = cat(1, posraw{idx_system}(:,1:2)', ones(1,size(posraw{idx_system},1)));
	    
        % Chain n iterations of the gait together, offsetting the start of each
        % by the displacement at the end of the previous gait, and avoiding a
        % double frame at the transition from one cycle to the next

        % Initialize the position data with an empty matrix
        posdata{idx_system} = [];

        % Iterate over the n gait cycles
        for i = 1:n_gaits

            % Skip first value if not the starting cycle (so that we don't get
            % a pause from a doubled frame at the changeover from one cycle to
            % the next
            if i == 1
                first_index = 1;
            else
                first_index = 2;
            end

            % Offset the displacements, using low-level implementation of SE(2)
            % action 

            % xy components (uses the start-pos to offset the starting xy
            % position so that it passes through the center of the frame at the
            % halfway point of the movie
            posdata_xy_augmented = ...
                powerm_pade(cyclic_displacement_m,...  % Take the displacement over one gait cycle
                    +start_pos...                      % Raise it to the negative power of the starting offset (so that we start back from the origin
                    +(i-1))...                         % For each subsequent gait, shift if forward by one gait cycle
                *posraw_xy_augmented(:,first_index:end); % Apply the starting offset for this cycle (calculated in the lines above) to the within-cycle displacements

            % theta component (always starts at zero orientation)		
            posdata_theta = ...
                posraw{idx_system}(first_index:end,3)'... % Take the theta values from the within-gait displacement
                +cyclic_displacement(3)*(i-1); % Add as many net theta changes as there are previous gait cycles

            % merge xy and theta data (stripping off ones from xy position at the same time)
            posdata_new = cat(1,posdata_xy_augmented(1:2,:),posdata_theta);

            % Concatenate the displacements
            posdata{idx_system} = cat(2,posdata{idx_system},posdata_new);
            
        end

			
        %return position data to column format
        posdata{idx_system} = posdata{idx_system}';
        
    end

end