function varargout = nicebox(varargin)
% Nicebox was originally a way of getting top edges on a plot without a
% full wireframe box. This behavior is now standard in MATLAB, so this code
% is deprecated and now just acts as a wrapper for BOX.

box(varargin{:})

end


% %NICEBOX draws an outline box on a set of axes, similar to the standard BOX
% %command, except that on a 3d plot, the three edges connecting to the point
% %closest to the camera will not be drawn. NICEBOX will suppress the lines
% %drawn by BOX, but restore them if toggled off
% %
% %NICEBOX('on') draws the nicebox lines for the current axis
% %NICEBOX('off') removes the nicebox lines for the current axis
% %NICEBOX or NICEBOX('') toggles the nicebox state of the current axes
% %NICEBOX('redraw') redraws the nicebox
% %NICEBOX(AX,...) uses AX instead of the current axes
% %
% %Lines drawn by nicebox will automatically update in response to
% %mouse-controlled 3d rotation. After other changes in the camera view
% %(such as with the VIEW command) or changes to the axis limits,
% %NICEBOX('redraw') should be used to update the lines.
% 
% 	%%%%%%%%%%%%
% 	%Handle the inputs
% 
% 	%no inputs specified
% 	switch nargin
% 
% 		case 0
% 
% 			%work on the current axis
% 			ax = gca;
% 
% 			%set the mode as 'toggle'
% 			mode = 'toggle';
% 
% 		case 1
% 
% 			%First check if the argument is a handle to a set of axes
% 			if any(ishandle(arg1)) && any(strcmpi('axes',get(arg1,'type')))
% 
% 				%set the target axis to arg1
% 				ax = arg1;
% 
% 				%set the mode to 'toggle'
% 				mode = 'toggle';
% 
% 			%Otherwise, check if the argument is a supported string
% 			elseif ischar(arg1) && any(strcmpi(arg1,{'on','off','','redraw'}))
% 
% 				%set the axis to the current one
% 				ax = gca;
% 
% 				%set the mode
% 				mode = arg1;
% 
% 				%if the mode is blank, change it to 'toggle'
% 				if isempty(mode)
% 					mode = 'toggle';
% 				end
% 
% 			else
% 
% 				error('Single argument to nicebox is not an axis handle or supported string')
% 
% 			end
% 
% 		case 2
% 
% 			%Check if the first argument is a handle to a set of axes
% 			if ishandle(arg1) && strcmpi('axes',get(arg1,'type'))
% 
% 				%set the target axis to arg1
% 				ax = arg1;
% 
% 			else
% 
% 				error('First argument to nicebox should be an axis handle')
% 
% 			end
% 
% 			%Check if the second argument is a supported string
% 			if ischar(arg2) && any(strcmpi(arg2,{'on','off','','redraw'}))
% 
% 				%set the mode
% 				mode = arg2;
% 
% 				%if the mode is blank, change it to 'toggle'
% 				if isempty(mode)
% 					mode = 'toggle';
% 				end
% 
% 			else
% 
% 				error('Second argument to nicebox should be a supported string')
% 
% 			end
% 
% 	end
% 
% 
% 	%%%%%%%%%%%
% 	%Get the nicebox state from the current axes
% 	udata = get(ax,'UserData');
% 	if isfield(udata,'nicedata')
% 		nicedata = udata.nicedata;
% 	else %if nicedata doesn't exist, start building it
% 		nicedata.onoff = 'off';
% 		nicedata.h = zeros(12,1);
% 		for i = 1:length(nicedata.h)
% 			nicedata.h(i) = ...
% 				line('Parent',ax,'XData',[],'YData',[],'ZData',[],'handleVisibility','off');
% 		end
% 		nicedata.rotate3d = rotate3d(ax);
% 		set(nicedata.rotate3d,'ActionPostCallback',@rerun_nicebox);
% 	end
% 
% 	%%%%%%%%%%%
% 	%Get the all the vertex points and make the connectivity grid;
% 
% 	Xlim = get(ax,'Xlim');
% 	Ylim = get(ax,'Ylim');
% 	Zlim = get(ax,'Zlim');
% 
% 	vertices = [...
% 		min(Xlim) min(Ylim) min(Zlim); %xyz
% 		max(Xlim) min(Ylim) min(Zlim); %Xyz
% 		min(Xlim) max(Ylim) min(Zlim); %xYz
% 		max(Xlim) max(Ylim) min(Zlim); %XYz
% 		min(Xlim) min(Ylim) max(Zlim); %xyZ
% 		max(Xlim) min(Ylim) max(Zlim); %XyZ
% 		min(Xlim) max(Ylim) max(Zlim); %xYZ
% 		max(Xlim) max(Ylim) max(Zlim)]; %XYZ
% 
% 	connectivity = [...
% 		1	2;	%min z plane
% 		1	3;
% 		2	4;
% 		3	4;
% 		1	5;	%vertical lines
% 		2	6;
% 		3	7;
% 		4	8;
% 		5	6;	%max Z plane
% 		5	7;
% 		6	8;
% 		7	8];
% 
% 	%%%%%%%%%
% 	%Get the point closest to the camera
% 
% 	%position of the camera
% 	campos = repmat(get(ax,'CameraPosition'),[8,1]);
% 
% 	%squared distance of the camera from vertex points
% 	camdist = sum((campos-vertices).^2,2);
% 
% 	%Get the index of the minimum-distance vertex
% 	[junk, camI] = min(camdist);
% 
% 
% 	%%%%%%%%%%
% 	%Make the lines that don't include the vertex closest to the camera
% 
% 	%remove lines including the camera vertex from the connectivity graph
% 	connectivity(((connectivity(:,1) == camI) | (connectivity(:,2) == camI)),:) = [];
% 
% 	%%
% 	%Turn the connectivity graph into an XYZ structue, with NaN values between
% 	%edges
% 
% 	%prime edge data structure
% 	edgedata = NaN(size(connectivity,1)*3,3);
% 	
% 	% fill in the edges
% 	for i = 1:size(connectivity,1)
% 		
% 		set(nicedata.h(i),'XData',vertices(connectivity(i,:),1)...
% 		,'YData',vertices(connectivity(i,:),2)...
% 		,'ZData',vertices(connectivity(i,:),3)...
% 		,'LineWidth',2*get(ax,'LineWidth'),'Color','k');
% 	
% 	end
% 		
% 
% % 	%copy in the vertex locations
% % 	for i = 1:size(connectivity,1)
% % 
% % 		edgedata(1+3*(i-1),:) = vertices(connectivity(i,1),:);
% % 		edgedata(2+3*(i-1),:) = vertices(connectivity(i,2),:);
% % 
% % 	end
% % 
% % 	%%%%%%%%
% % 	%Set the line data for the box
% % 
% % 	set(nicedata.h,'XData',edgedata(:,1)...
% % 		,'YData',edgedata(:,2)...
% % 		,'ZData',edgedata(:,3)...
% % 		,'LineWidth',2*get(ax,'LineWidth'),'Color','k')
% 
% 
% 
% 
% 	%%%%%%%
% 	%Set the visibility of the nicebox
% 
% 	%first, include the effect of any mode changes
% 	switch mode
% 
% 		case 'on'
% 
% 			nicedata.onoff = 'on';
% 
% 		case 'off'
% 
% 			nicedata.onoff = 'off';
% 
% 		case 'toggle'
% 
% 			switch nicedata.onoff
% 
% 				case 'on'
% 
% 					nicedata.onoff = 'off';
% 
% 				case 'off'
% 
% 					nicedata.onoff = 'on';
% 
% 			end
% 
% 		case 'redraw'
% 
% 			%do not change the onoff mode
% 
% 		otherwise
% 
% 			error('Unknown mode string')
% 
% 	end
% 
% 
% 	%Now use the onoff value to set the visibility
% 	set(nicedata.h,'Visible',nicedata.onoff)
% 
% 
% 
% 	%%%%%%%%
% 	%put the nicedata into the userdata of the axes
% 	udata.nicedata = nicedata;
% 	set(ax,'UserData',udata)
% 	
% 	%return the nicebox state if asked
% 	if nargout == 1;
% 		varargout = {nicedata.onoff};
% 	else
% 		varargout = {};
% 	end
% 
% end
% 
% 
% %Callback function to redraw the outline when the axes are manually rotated
% function rerun_nicebox(obj,evd) %#ok<INUSL>
% 
% 	%Call nicebox on the associated axes
% 	nicebox(evd.Axes,'redraw')
% 
% end