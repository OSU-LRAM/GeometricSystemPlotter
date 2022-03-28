function animate_locomotor(export,info_needed)
    
    % Look up the geometry specification for this system:
    sysfile = fullfile(info_needed.datapath, ['sysf_', info_needed.current_system2, '_calc.mat']);
    load(sysfile,'s')
	info_needed.s = s;
    
  	% Declare a directory name file names for the movies
	destination_root = fullfile(info_needed.datapath,'Animation',[info_needed.current_system2, '__' info_needed.current_shch2]);
	destination_list = fieldnames(info_needed.Movies);
            
        if strcmp(info_needed.Coordinates,'minperturbation_coords')
            destinationsuffix = 'min_perturb';
        else
            destinationsuffix = 'original_coords';
        end

	destination = cellfun(@(x) fullfile(destination_root,[info_needed.current_system2, '__' info_needed.current_shch2 '__' x '__' destinationsuffix]),destination_list,'UniformOutput',false); 

    % Flag movies that should be generated, don't skip any movies
    export_list = cellfun(@(x) info_needed.Movies.(x), destination_list);
    skip_list = zeros(size(export_list));
    info_needed.export_list = export_list;
    
    %%%
    % Now cull the destination list
    destination(~export_list) = [];
    
    %%%%%%%%
	% Create animation elements, and store them in the frame_info structure
	frame_info{1} = create_elements(info_needed); %setup function, defined below in file
	frame_info = [frame_info;create_elements_shapespace(info_needed)]; % set up shapespace plots

	
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
    frame_info{1}.s = s;

	% Animate the movie
	[frame_info, endframe]...
		= sysplotter_animation(frame_gen_function,frame_info,timing,destination,export_list,skip_list,0);

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
        'XLim',[-1,1]*.7/.45*info_needed.s.geometry.length,...  % Axes scaled to system scale
        'YLim',[-1,1]*.7*info_needed.s.geometry.length);
 	set(h.ax,'Position',[0 0 1 1])                  % Make the axes fill the whole window

    data_source = info_needed.datapath;
	system_name = info_needed.current_system2;
	gait_name = info_needed.current_shch2;

	%Create the locomotor
	h.robot = create_locomotor(...
        h.ax,...
        data_source,...
        system_name,...
        info_needed);
	
	% Extract the position and shape data
	gaitfile = fullfile(data_source,['sysf_' system_name '__shchf_' gait_name]);
    load(gaitfile);
	shaperaw = p.phi_locus_full{1}.shape;
    
    if strcmp(info_needed.Coordinates,'minperturbation_coords')
        posraw = p.G_locus_full{1}.G_opt;
        h.drawing_baseframe = info_needed.current_system2;
    else
        posraw = p.G_locus_full{1}.G;
    end

	
	
	% Specify the number of times to repeat the cycle, and the number of
	% negative cycle-displacements to start from
	n_gaits = info_needed.Number_gaits;
	start_pos = -(n_gaits/2);
	
	% get full configuration history
	[h.shapedata, h.posdata] = gait_concatenator(shaperaw,posraw,n_gaits,start_pos);
	
        load('sysplotter_config','Colorset')

    
	%%%%%%%%
	% Create the tracer objects
	h.tld = line('Parent',h.ax,'XData',[],'YData',[],'ZData',[],'LineWidth',5,'Color',Colorset.spot,'Marker','o','MarkerFaceColor',Colorset.spot,'MarkerSize',5,'LineStyle','none');
	h.tl = line('Parent',h.ax,'XData',[],'YData',[],'ZData',[],'LineWidth',5,'Color',Colorset.spot); %Tracer line for translation

	
end

function h = create_elements_shapespace(info_needed)

    % Pull information out of info_needed
    s = info_needed.s;
    export_list = info_needed.export_list;

	% Create lists to iterate along for shape space axis creation, and then
	% trim out movies that aren't being made
	fignums = [171; 172; 173; 181; 182; 183];
	fignums(~export_list(2:end)) = [];
    
	% Name the plots to draw into those windows, and trim them to just the
	% plots being made
	fignames = {'CCF X';'CCF Y';'CCF Theta'; 'Avec X';'Avec Y';'Avec Theta'};
    
        if strcmp(info_needed.Coordinates,'minperturbation_coords')
            fignamessuffix = ' min. perturbation coords';
        else
            fignamessuffix = ' original coords';
        end
        fignames = cellfun(@(x) [x fignamessuffix],fignames,'UniformOutput',false);

    fignames(~export_list(2:end)) = [];

    %%%%%
    % Internal names used by sys_draw
    figinternal = {'X';'Y';'T';'X';'Y';'T'};
        if strcmp(info_needed.Coordinates,'minperturbation_coords')
            figinternalsuffix = 'opt';
        else
            figinternalsuffix = [];
        end
        figinternal = cellfun(@(x) {[x figinternalsuffix]},figinternal,'UniformOutput',false);
        
    figinternal(~export_list(2:end)) = [];
     
    %%%%%%
    % Category of plots to make
    figcategory = {'CCF';'CCF';'CCF';'vfield';'vfield';'vfield'};
    figcategory(~export_list(2:end)) = [];
	
    %%%%%%
	% Generate axes with consistent figure numbers
	ax = zeros(size(fignums));
    h = cell(size(fignums));
    underlying_function = h;
	for idx = 1:length(fignums)
        
            f = figure(fignums(idx));
            clf(f);
            ax(idx) = axes('Parent',fignums(idx));
            set(f,...
                'renderer','painters',...
                'name',[info_needed.current_system2,' ',fignames{idx}])

            h{idx,1}.f = f; 
            h{idx}.ax = ax(idx);
            
            
            % Identify the function with the drawing code
            underlying_function{idx} = ...
                fullfile(info_needed.sysplotterpath,'sys_draw_fcns',[figcategory{idx} '_draw']);        
		
	end
	

	% Build the structure that sysplotter uses internally to call for plots
	plot_info = struct(...
        'axes',num2cell(ax),...
        'components',figinternal,...
		'category',figcategory,...
        'stretch',0,...
        'stretchpath', 'null',...
		'CCFtype','DA',...
        'stretch_name','none',...
        'style','contour');

	% Identify the data file to draw the plot information from
	data_source = info_needed.datapath;
	system_name = info_needed.current_system2;
	load(fullfile(data_source,['sysf_' system_name '_calc']),'s');
	



    % Call the drawing function

    resolution.vector = s.density.vector;

    resolution.scalar = s.density.scalar;
    
    % TODO: this shouldn't be hardcoded...
    resolution.vector_range = s.grid_range;

    resolution.scalar_range = s.grid_range;

    % Call the drawing function
    for idx = 1:numel(plot_info)
        hh(idx,1) = absolute_feval(underlying_function{idx},s,[],plot_info(idx),[],'null',resolution);
    end
	
	% Remove click callbacks on axes created
    if exist('hh','var')
        for i = 1:numel(hh)
            % Remove the callbacks on clicks
            set(hh(i).axes,'ButtonDownFcn',[])

            % Label axes
            set(get(hh(i).axes,'xlabel'),'string','$a_{1}$')
            set(get(hh(i).axes,'ylabel'),'string','$a_{2}$')
        end
    end
    
    load('sysplotter_config','Colorset')
    
    % Now create a tracer-and-dot on each axis
	for i = 1:numel(h)
		h{i}.tl = line('Parent',h{i}.ax,'XData',[],'YData',[],'ZData',[],'LineWidth',7,'Color',Colorset.spot); %Tracer line
		h{i}.tld = line('Parent',h{i}.ax,'XData',[],'YData',[],'ZData',[],'LineWidth',5,'Color',Colorset.spot,'Marker','o','MarkerFaceColor',Colorset.spot,'MarkerSize',10,'LineStyle','none');
	end

end


%%%%%%%%%%%%%%%%%%
% Frame content to draw during the gait
function frame_info = execute_gait(frame_info,tau)

	% Call out the frame_info pieces that contain robot verse shapespace
	% data
	irobot = 1;
	ishape = 2:size(frame_info);

	% Generate a core timing vector from zero to 1, with as many entries as
	% data points in the kinematic history
	timing_base = linspace(0,1,size(frame_info{irobot}.shapedata,1));
	
	%%%%%%
    % Get the system position and shape at fractional time tau
    
    % Interpolate the shape variables at the specified time
    config.shape = zeros(size(frame_info{irobot}.shapedata,2),1);
    for idx = 1:numel(config.shape)
        config.shape(idx) = interp1(timing_base,frame_info{irobot}.shapedata(:,idx),tau);
    end
    
    % Interpolate the position variables at the specified time
	config.x = interp1(timing_base,frame_info{irobot}.posdata(:,1),tau);
	config.y = interp1(timing_base,frame_info{irobot}.posdata(:,2),tau);
	config.theta = interp1(timing_base,frame_info{irobot}.posdata(:,3),tau);
	
    %%%%
	% Place the locomotor at the position
    
    % If a drawing baseframe has been specified, append the original
    % system baseframe with the one that has been specified
    if isfield(frame_info{irobot},'drawing_baseframe') && ~isfield(frame_info{irobot},'drawing_baseframe_inserted')
        if ~iscell(frame_info{irobot}.s.geometry.baseframe)
            frame_info{irobot}.s.geometry.baseframe ...
                = {frame_info{irobot}.s.geometry.baseframe};
        end
        if ~iscell(frame_info{irobot}.s.geometry.baseframe)
            frame_info{irobot}.drawing_baseframe ...
                = {frame_info{irobot}.drawing_baseframe};
        end
        frame_info{irobot}.s.geometry.baseframe = [frame_info{irobot}.s.geometry.baseframe frame_info{irobot}.drawing_baseframe];
        frame_info{irobot}.drawing_baseframe_inserted = 1;
    end
    
    % Use the configuration to place the locomotor
  	frame_info{irobot}.robot = place_locomotor(frame_info{irobot}.robot,config,frame_info{irobot}.s);
	
	% draw the locomotor
	frame_info{irobot}.robot = draw_locomotor(frame_info{irobot}.robot,0);
	
	% Draw the tracer dot
	set(frame_info{irobot}.tld,'XData',[frame_info{irobot}.posdata(1,1) config.x],'YData',[frame_info{irobot}.posdata(1,2) config.y],'ZData',[5,5]);
	
	% Draw the tracer
	set(frame_info{irobot}.tl,'XData',frame_info{irobot}.posdata(1:floor(tau*end),1),'YData',frame_info{irobot}.posdata(1:floor(tau*end),2),'ZData',5*ones(size(frame_info{irobot}.posdata(1:floor(tau*end),1))));

	% Declare a print method (in this case, print 150dpi png files of
	frame_info{irobot}.printmethod = @(dest) print(frame_info{irobot}.f,'-r 150','-painters',dest);

	%%%%%%%%%%%%%%%%
	% Shapespace drawing stuff
	
    % If there are shape axes, draw the curve from the start time up to the
    % current time onto the axes, and put dots at the beginning and end
	shapeh = [frame_info{ishape}];
    if ~isempty(shapeh)
        set([shapeh.tl],'XData',frame_info{irobot}.shapedata(1:ceil(tau*end),1),'YData',frame_info{irobot}.shapedata(1:ceil(tau*end),2),'ZData',5*ones(size(frame_info{irobot}.shapedata(1:ceil(tau*end),1))));
        set([shapeh.tld],'XData',[frame_info{irobot}.shapedata(1,1) config.shape(1)],'YData',[frame_info{irobot}.shapedata(1,2) config.shape(2)],'ZData',[5,5]);
    end
    
    % Set the print method for the frames
	for i = 1:numel(ishape)
		frame_info{ishape(i)}.printmethod = @(dest) print(frame_info{ishape(i)}.f,dest);
	end
	

end

%%%%%%
function [shapedata, posdata] = gait_concatenator(shaperaw,posraw,n_gaits,start_pos)


	
	% Store the shape and data to the frameinfo structure
	shapedata = cat(1,shaperaw(1,:),repmat(shaperaw(2:end,:),[n_gaits,1]));
	
	%%%%%%
	% Concatenate the gait displacements together.
	
	% Get the displacement over one cycle
	cyclic_displacement = posraw(end,:);
	
	% convert this cyclic displacement into an SE(2) matrix
	cyclic_displacement_m = vec_to_mat_SE2(cyclic_displacement);
	
    % Add a set of ones onto the xy values, so that they can be
    % translated/rotated by SE(2) matrices
	posraw_xy_augmented = cat(1, posraw(:,1:2)', ones(1,size(posraw,1)));
	    
	% Chain n iterations of the gait together, offsetting the start of each
	% by the displacement at the end of the previous gait, and avoiding a
	% double frame at the transition from one cycle to the next
    
    % Initialize the position data with an empty matrix
	posdata = [];

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
            posraw(first_index:end,3)'... % Take the theta values from the within-gait displacement
            +cyclic_displacement(3)*(i-1); % Add as many net theta changes as there are previous gait cycles
		
		% merge xy and theta data (stripping off ones from xy position at the same time)
		posdata_new = cat(1,posdata_xy_augmented(1:2,:),posdata_theta);
		
		% Concatenate the displacements
		posdata = cat(2,posdata,posdata_new);        

		
	end
	
	%return position data to column format
	posdata = posdata';

end

%%%%%%
function label_shapespace_axes(ax,z_text)
%put the x and y labels on the plots. optional argument 'z_text' sets a
%label for the z axis

    %Default text labels
    x_text = '$\alpha_1$';
    y_text = '$\alpha_2$';
    
    %Formatting options
    format_list = {'FontName','Times','Interpreter','latex','FontSize',30};
    
    %apply the label
    xlabel(ax,x_text,format_list{:});
    ylabel(ax,y_text,format_list{:});
    
    %if there's a third dimension, label it
    if exist('z_text','var')
        
        zlabel(ax,z_text,format_list{:});
        
    end
        
    
end



%%%%%%
function set_tics_shapespace(ax,s)
%place the tic marks
    
    %set the tic fontsize
    set(ax,'FontSize',20,'FontName','GillSans')
    set(ax,'XTick',s.tic_locs.x,'YTick',s.tic_locs.y)

end