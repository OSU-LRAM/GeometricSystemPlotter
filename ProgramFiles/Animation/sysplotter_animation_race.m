function [frame_info, endframe] = sysplotter_animation_race(frame_gen_function,frame_info,timing,destination,export,skip,startframe)
% An animation framework for writing a series of frames into a directory,
% at which point ffmpeg, QuicktimePro, or the like can be used to turn them
% into a movie.
%
%
% ================
%
% Inputs are:
%
% frame_gen_function: A function that takes the frame_info structure and
% normalized time (from 0 to 1) and draws
% draws the animation frame for this movie at that time, returning a frame_info
% structure that may be updated with any new information about the frame
% content. *Note that frame_gen_function should indicate the desired printing
% mode for the frame in the frame_info structure as noted below*
%
% frame_info: A structure containing:
%	frame_info.printmethod = the print function, with all parameters set
%	 except for destination, which will be set based on the 'destination'
%	 parameter here
%	any other fields (such as axes, or information about the content of the
%	 previous frame) used within frame_gen_function.
%
% timing: A structure containing the timing information for the animation:
%	timing.duration: How many seconds the animation should run for
%	timing.fps: How many frames per second in the animation
%	timing.pacing: Handle to a pacing function. This should be something
%	like 
%      @(y) linspace(0,1,y)
%   i.e. taking a number of frames to plot, and sampling them across a
%   specified number of frames according to a desired distribution. The
%   interval will be normalized to a unit interval post-sampling.
% 
% destination: string indicating the path (local or global) of the
% directory in which to place the images.
%
% export: optional argument, boolean for whether to export the frames or
% just draw them. Defaults to 0 (no export)
%
% skip: optional argument, will draw only the last frame of the animation.
% defaults to 0 (do not skip the animation)
%
% startframe: optional argument to set the starting frame. Defaults to 1.
%
% ========================
%
% Outputs are:
%
% frame_info: last frame_info output by the frame_gen_function
%
% endframe: last frame number written by the program
%
%
% ========================
% Copyright Ross L. Hatton, 2011

    % Create the object for saving the file
%     writerObj = VideoWriter('swimmer.avi');


	%%%%
	% Handle default arguments


	% Set the export to default if necessary
	if ~exist('export','var')
		export = 0;
	end
	
	% Set the skip to default if necessary
	if ~exist('skip','var')
		skip = 0;
	end

	
	% Set the startframe to default if necessary
	if ~exist('startframe','var')
		startframe = 1;
    end
    
%     % view will be used like export, but for which animations you want to
%     % have play live (unexported, in a matlab figure)
%     view = export; % For now, all exported animations are also played after.
%     % animations that are either viewed or exported need frames calculated
%     calc = arrayfun(@or,export,view);
%     % warning: export and view must be the same size

	
	% If the destination is a single string, then wrap it in a cell array
	if isstr(destination)
		destination = {destination};
	end
	
    % Set output video filename
    if isfield(frame_info{1}, 'sysname') && isfield(frame_info{1}, 'pathname')
        outputfname = [frame_info{1}.sysname '_' frame_info{1}.pathname];
    else
        outputfname = 'gait_animation';
    end
    
	% Globalize the path if necessary (some printing commands may require
	% global destinations). Checks to see if destination starts with an
	% absolute-path specifier, and prepends the working directory if this
	% is not present
    Animation_dir = 'Animation';
    
	if ~ ( (ispc && strcmp(':\',destination{1}(2:3))) || (isunix && strncmp('/',destination{1},1) ) )
		
		destination = cellfun(@(x)fullfile(pwd, x, Animation_dir),destination,'UniformOutput',false);
			
	end
	
	% Only ensure directories are present if export is set
%	for i = 1:numel(export)
    for i = 1:numel(destination)
		
%		if export(i)
			
			% Ensure that destination directory exists
            [parent,movie] = fileparts(destination{i});
			if ~exist(parent,'dir')
				mkdir(parent)
			end


		
%		end
		
    end

	
	%%%%%%%%%%%%%%%%%
	% Work out the timing. Get the framepoints (fractions within the span of
	% the movie, normalized to a unit interval)
	
	% Generate a set of points from the timing function
	framepoints_raw = timing.pacing(round(timing.duration * timing.fps));
	
	% Put these points on a zero-to-one scale (necessary if, e.g., the
	% timing vector is logarithmic)
	framepoints_zerostart = framepoints_raw-framepoints_raw(1);
	framepoints = framepoints_zerostart/framepoints_zerostart(end);
	
	%%%%%%%%%%%%%%
	% Write the frames
	
	% Prime the framenumber
	framecounter = startframe;
	
    % Prepare datastructure to store frames
    frame_info = frame_gen_function(frame_info,framepoints(1));

    firstframe = frame_info{1}.printmethod('-RGBImage');
    for j = 1:numel(destination)
        %if calc(j)
            F{j} = struct('cdata',num2cell(zeros([length(framepoints) size(firstframe)],class(firstframe)),[2,3,4]),'colormap',[]);
        %end
    end
        
	% Print out the frames
	if ~skip							% If skipping this movie, don't draw frames
		for i = 1:length(framepoints)

			% Draw the frame, and get the information about how to print the
			% frame
			frame_info = frame_gen_function(frame_info,framepoints(i));
			
			% Ensure that all elements of the figure are updated.
			drawnow

			% Turn framecounter into a string. String is set up to prepend
			% zeros to turn this into a 8-digit number.
% 			frame_str = num2str(framecounter);
% 			frame_str = ['fr' repmat('0',1,8-length(frame_str)) frame_str]; %#ok<AGROW>

% 			% Print out the frame into the newframes directory if the
% 			% export flag is set
            % Save animation datastructures for all 
			for j = 1%:numel(destination)
				
%				if calc(j)
            
                    % Save current frame
                    % Compatability with old code warning: Remove the image type
                    % input from the printmethod function and it should work.
%                     F{j}(i).cdata = rgb2ind(frame_info{j}.printmethod('-RGBImage'),colormap,'nodither');
                    F{j}(i).cdata = frame_info{j}.printmethod('-RGBImage');
                    % TODO: change F into a cell array of these sorts
                    % of things so that the gait tracing animations can be
                    % saved in it too
                    
% 					destination_str = fullfile(destination{j},'newframes',frame_str);
% 					frame_info{1}.printmethod(destination_str);
					
%				end
			end

			% Update the frame counter
			framecounter = startframe+i;
            
%             % Open the object
%             open(writerObj);
%             % Get the frame
%             M = getframe(171);
%             % Make the AVi file
%             writeVideo(writerObj,M)

		end
		
	else % execute the last frame, to update the image, and adjust the framecounter
		
		% Draw the last frame of the movie
		frame_info = frame_gen_function(frame_info,framepoints(end));
		
		% clear the drawing queue
		drawnow
		
		% Ensure that the framecounter is the same as if the movie hadn't
		% been skipped
		framecounter = startframe+length(framepoints);
	end
	
	% If the export flag was set, remove the old frames, and move the new
	% frames to the main directory
% 	for i = 1:numel(export)
% 		if export(i)
% 			try
% 				rmdir(fullfile(destination{i},'oldframes'));
% 			end
% 
% 			movefile(fullfile(destination{i},'newframes','fr*'),destination{i});
% 
% 		end
% 	end
	
	% Return the last frame printed
	endframe = framecounter-1;

    
    warning('off','MATLAB:audiovideo:VideoWriter:mp4FramePadded')
    % Export:
    for j = 1%:numel(destination)
        %if export(j)

            if ispc || ismac
                videoformat = 'MPEG-4';
            else
                videoformat = 'Motion JPEG AVI';
            end

            v = VideoWriter(destination{j},videoformat);
            v.FrameRate = timing.fps;
            open(v)
            writeVideo(v,F{j})
            close(v)
                
        %end
    end

%     % Close the figures used during animation making:
%     fignums = [17; 171; 172; 173]; % Hardcoded from animate_backbone...
%     for fig = fignums
%        close(fig) 
%     end
    
%     % Play:
%     for j = 1:numel(calc)
%         if view(j)
%             % Create a figure and resize to fit the animation
%             animfig=figure(17+j);
%             clf(animfig)
%             animsize=size(F{j}(end).cdata);
%             animfig.Position(3:4)=[animsize(2) animsize(1)];
%             movegui(animfig) % ensure after resizing that figure is still onscreen	
%     % Play the movie:
% %     Maybe resize in order to show live animation smaller?
% %     Can't get the code to work right now...
% %     Fsm{j}.cdata = arrayfun(@(x) imresize(x.cdata,0.5),F{j},'UniformOutput',false);
% %     Fsm{j}.colormap = F{j}.colormap;
% %             movie(animfig,Fsm{j},1,timing.fps)
%             movie(animfig,F{j},1,timing.fps)
%         end
%     end
	
end