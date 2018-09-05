function animate_4_links_swimmer(export,info_needed)

	%%%%%%%%
    addpath('robot_drawing_tools/robot_drawing_tools')
	% Create animation elements, and store them in the frame_info structure
	frame_info{1} = create_elements(info_needed); %setup function, defined below in file
	frame_info = [frame_info;create_elements_shapespace(info_needed)]; % set up shapespace plots

	
	% Designate animation function
	frame_gen_function...
		= @execute_gait; % frame function, defined below in file

	% Declare timing
	timing.duration = info_needed.Timeduration; % three second animation
	timing.fps = info_needed.Framerate;     % create frames for 15 fps animation
	timing.pacing = @(y) linspace(0,1,y); % Use a soft start and end, using the included softstart function

	% Declare a directory name in which to place files
	destination_root = info_needed.current_system2;
	destination_list = {'snake','x','y','theta'};
	destination = cellfun(@(x) fullfile(destination_root,x),destination_list,'UniformOutput',false); 

	% Animate the movie
	[frame_info, endframe]...
		= animation(frame_gen_function,frame_info,timing,destination,[0,0],[0,0],0);

end

function h = create_elements(info_needed)

	h.f = figure(17);                            % Designate a figure for this animation
	clf(h.f)                                     % Clear this figure
	set(h.f,'position',[600 1500 800 350],'paperposition',[0 0 10 4.5],'papersize',[10 4.5],'Color','w');

	h.ax = axes('Parent',h.f);                   % Create axes for the plot
	axis(h.ax,'equal','off');
	set(h.ax,'XLim',[-5,5]*.8/2.5,'YLim',[-2.25,2.25]*.8/2.5);
	set(h.ax,'Position',[0 0 1 1])


% 	data_source = '/Users/rlhatton/Documents/MATLAB/GeometricSystemPlotter/UserFiles/v4/rlhatton/sysplotter_data/v4.0/';
    
    data_source = info_needed.UserFile_path;
	system_name = info_needed.current_system2;
	gait_name = info_needed.current_shch2;

	%Create the robot
	h.robot = create_4_links_swimmer(h.ax,data_source,system_name);
	
	% Extract the position and shape data
	load(fullfile(data_source,['sysf_' system_name '__shchf_' gait_name]));
	shaperaw = p.phi_locus_full{1}.shape;
	posraw = p.G_locus_full{1}.G_opt;
	
	% Specify the number of times to repeat the cycle, and the number of
	% negative cycle-displacements to start from
	n_gaits = info_needed.Number_gaits;
	start_pos = -2;%(n_gaits/2)-1;
	
	% get full configuration history
	[h.shapedata, h.posdata] = gait_concatenator(shaperaw,posraw,n_gaits,start_pos);
	
	%%%%%%%%
	% Create the tracer objects
	h.tld = line('Parent',h.ax,'XData',[],'YData',[],'ZData',[],'LineWidth',5,'Color',[234 14 30]/255,'Marker','o','MarkerFaceColor',[234 14 30]/255,'MarkerSize',5,'LineStyle','none');
	h.tl = line('Parent',h.ax,'XData',[],'YData',[],'ZData',[],'LineWidth',5,'Color',[234 14 30]/255); %Tracer line for translation

	

end

function h = create_elements_shapespace(info_needed)

    addpath('\Users\hosse\Documents\MATLAB\Dr Hatton Snake Robot\GeometricSystemPlotter\ProgramFiles\v4.1\Utilities')


	% Create lists to iterate along for shape space axis creation
	fignums = [171; 172; 173];
	
	
	% Generate axes with consistent figure numbers
	ax = zeros(size(fignums));
	for i = 1:length(fignums)
		f = figure(fignums(i));
		clf(f);
		ax(i) = axes('Parent',fignums(i));
		
		h{i,1}.f = f; %#ok<AGROW>
		h{i}.ax = ax(i); %#ok<AGROW>
		
	end
	
	% Name the plots I want to draw into those windows
	figplots = {'Xopt','Yopt','Topt'};

	% Build the structure that sysplotter uses internally to call for plots
	plot_info = struct('axes',{ax},'components',{figplots},...
		'category',{'vfield'},'stretch',{0},'stretchpath', {'null'},...
		'CCFtype',{'DA'},'style',{'contour'});

	% Identify the data file to draw the plot information from
	data_source = info_needed.UserFile_path;
	system_name = info_needed.current_system2;
	load(fullfile(data_source,['sysf_' system_name '_calc']));
	
	% Identify the function with the drawing code
	vfieldfunction = fullfile(info_needed.current_directory,'\sys_draw_fcns\CCF_draw');

    % Call the drawing function

    resolution.vector = [11 11];

    resolution.scalar = [21 21];

    resolution.vector_range = [-2.5000 2.5000 -2.5000 2.5000];

    resolution.scalar_range = [-2.5000 2.5000 -2.5000 2.5000];

%     hh = absolute_feval(vfieldfunction,s,[],plot_info,[],'null',[],resolution);
% 	
% 	% Remove the callbacks on clicks
% 	set(hh.axes,'ButtonDownFcn',{})
	
	% Now create a tracer-and-dot on each axis
	for i = 1:numel(h)
		h{i}.tl = line('Parent',h{i}.ax,'XData',[],'YData',[],'ZData',[],'LineWidth',7,'Color',[234 14 30]/255); %Tracer line
		h{i}.tld = line('Parent',h{i}.ax,'XData',[],'YData',[],'ZData',[],'LineWidth',5,'Color',[234 14 30]/255,'Marker','o','MarkerFaceColor',[234 14 30]/255,'MarkerSize',10,'LineStyle','none');
	end

end


%%%%%%%%%%%%%%%%%%
% Frame content specification for dynamically drawing the cosine function,
% at time tau on a range of 0 to 1
function frame_info = execute_gait(frame_info,tau)

	% Call out the frame_info pieces that contain robot verse shapespace
	% data
	irobot = 1;
	ishape = 2:4;

	% Generate a core timing vector from zero to 1, with as many entries as
	% data points in the kinematic history
	timing_base = linspace(0,1,size(frame_info{irobot}.shapedata,1));
	
	% Get the system position and shape at fractional time tau
	config.alpha_1 = interp1(timing_base,frame_info{irobot}.shapedata(:,1),tau);
	config.alpha_2 = interp1(timing_base,frame_info{irobot}.shapedata(:,2),tau);
    config.alpha_3 = interp1(timing_base,frame_info{irobot}.shapedata(:,3),tau);
	
	config.x = interp1(timing_base,frame_info{irobot}.posdata(:,1),tau);
	config.y = interp1(timing_base,frame_info{irobot}.posdata(:,2),tau);
	config.theta = interp1(timing_base,frame_info{irobot}.posdata(:,3),tau);
	
	% Place the robot at the position
	frame_info{irobot}.robot = place_4_links_swimmer(frame_info{irobot}.robot,config,'mean');
	
	% draw the robot using the swimmer graphics
	frame_info{irobot}.robot = draw_swimmer(frame_info{irobot}.robot,1);
	
	% Draw the tracer dot
	set(frame_info{irobot}.tld,'XData',[frame_info{irobot}.posdata(1,1) config.x],'YData',[frame_info{irobot}.posdata(1,2) config.y],'ZData',[5,5]);
	
	% Draw the tracer
	set(frame_info{irobot}.tl,'XData',frame_info{irobot}.posdata(1:floor(tau*end),1),'YData',frame_info{irobot}.posdata(1:floor(tau*end),2),'ZData',5*ones(size(frame_info{irobot}.posdata(1:floor(tau*end),1))));

	% Declare a print method (in this case, print 150dpi png files of
	% frame_info{irobot}.f, using the Painters renderer)
	frame_info{irobot}.printmethod = @(dest) print(frame_info{irobot}.f,'-dpng','-r 150','-painters',dest);

	%%%%%%%%%%%%%%%%
	% Shapespace drawing stuff
	
	shapeh = [frame_info{ishape}];
	set([shapeh.tl],'XData',frame_info{irobot}.shapedata(1:ceil(tau*end),1),'YData',frame_info{irobot}.shapedata(1:ceil(tau*end),2),'ZData',5*ones(size(frame_info{irobot}.shapedata(1:ceil(tau*end),1))));
	set([shapeh.tld],'XData',[frame_info{irobot}.shapedata(1,1) config.alpha_1],'YData',[frame_info{irobot}.shapedata(1,2) config.alpha_2],'ZData',[5,5]);

	for i = 1:numel(ishape)
		frame_info{ishape(i)}.printmethod = @(dest) printeps_gillsans(frame_info{ishape(i)}.f,dest);
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
	cyclic_displacement_m = vec_to_mat_se2(cyclic_displacement);
	
	% the first set of displacements is the raw data, modified by the
	% cyclic displacement to the start_pos power. Each subsequent set is
	% similarly defined, but starting from the second value, and
	% incrementing the offset power by one cycle
	posraw_xy_augmented = cat(1, posraw(:,1:2)', ones(1,size(posraw,1)));
	
	% will build matrix incrementally
	posdata = [];
	for i = 1:n_gaits
		
		% Skip first value if not the starting cycle
		if i == 1;
			first_index = 1;
		else
			first_index = 2;
		end
		
		% Offset the displacements, using low-level implementation of SE(2)
		% action 
		
		% xy components
		posdata_xy_augmented = powerm_pade(cyclic_displacement_m,start_pos+(i-1))...
			*posraw_xy_augmented(:,first_index:end);
		
		% theta component		
		posdata_theta = posraw(first_index:end,3)'-posraw(end,3)*(start_pos+(i-1));
		
		% merge xy and theta data (stripping off augmentation at same time)
		posdata_new = cat(1,posdata_xy_augmented(1:2,:),posdata_theta);
		
		% Concatenate the displacements
		posdata = cat(2,posdata,posdata_new);
		
	end
	
	%return to column format
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
