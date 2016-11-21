function painter_layer(f,layers,destination,printoptions)
%Print out an image to a set of layers, then re-compose them. This gets
%around the layering problems with the Painters renderer. Each layer's
%objects should be a vector of handles, and all the layers should be
%grouped into a 1-d cell array, with background layers first

%make sure all my scripts are available
setenv('PATH','/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/usr/X11/bin:/usr/texbin:/opt/local/bin:/opt/local/sbin:/Users/rlhatton/scripts')

%make sure that the tempdir for building the composite image exists and is
%empty

recycle('on'); %move files to trash rather than rm-ing them

if ~exist('~/Temp/painted_layer/','dir')
	mkdir('~/Temp/painted_layer/');
else
	delete('~/Temp/painted_layer/*')
end


%%%%%%%%%%%%%%%%%%

%Record the nicebox state of any axes in the figure
for i = 1:length(layers)
		
	%verify that the elements for this layer are handles
	if ~all(ishandle(layers{i}))
		error(['There is a nonhandle in layer ' num2str(i) '.'])
	end
	
	%check the elements of the layer for axes
	axes_found = strcmp('axes',get(layers{i},'type'));
	for j = 1:length(axes_found)
		
		if axes_found(j)
			
			nicebox_state{i}{j} = nicebox(layers{i}(j),'redraw');
			
		end
		
	end
	
end

%Work down the specified layers
for i = 1:length(layers)
		
	%Make the elements for this layer visible, and all others invisible
	set(layers{i},'Visible','on')
	
	%bring nicebox to the appropriate state for any axes in the layer

	%check the elements of the layer for axes
	axes_found = strcmp('axes',get(layers{i},'type'));
	if any(axes_found)
		for j = 1:length(axes_found)

			if axes_found(j)

				nicebox(layers{i}(j),nicebox_state{i}{j});

			end

		end
	end
		
	
	%make all the other layers invisible
	j_list = 1:length(layers);
	j_list(i) = [];
	for j = j_list
				
		%Make the elements for this layer visible, and all others invisible
		set(layers{j},'Visible','off')
		
		%check the elements of the layer for axes
		axes_found_b = strcmp('axes',get(layers{j},'type'));
		if any(axes_found_b)
			for k = 1:length(axes_found_b)

				if axes_found_b(k)

					nicebox(layers{j}(k),'off');

				end

			end
		end
	

	end
	
	
	
	% Set the renderer
	%check the elements of the layer for surfaces
	surf_found = strcmp('surface',get(layers{i},'type'));
	if any(surf_found)
		renderer_string = {'-painters'};
	else
		renderer_string = {'-painters'};
	end
		
	
	%print this layer
	print(f,'-depsc2','-loose',renderer_string{:},printoptions{:},['~/Temp/painted_layer/' num2str(i)])		
	%epsfontsubs(f,['~/Temp/painted_layer/' num2str(i)],printoptions{:})
	
end

%Restore the state of the figure
%Work down the specified layers
for i = 1:length(layers)
		
	%Make the elements for this layer visible, and all others invisible
	set(layers{i},'Visible','on')
	
	%bring nicebox to the appropriate state for any axes in the layer

	%check the elements of the layer for axes
	axes_found = strcmp('axes',get(layers{i},'type'));
	if any(axes_found)
		for j = 1:length(axes_found)

			if axes_found(j)

				nicebox(layers{i}(j),nicebox_state{i}{j});

			end

		end
	end
end

%%%%%%%%%%%%%%%%%%%%
%compose the layers

%build up the arguments for epscompose
filenames = [num2str((1:length(layers))'),repmat('.eps 0 0 1 0 ',length(layers),1)]';

unix(['cd ~/Temp/painted_layer/; epscompose ' (filenames(:))' ' > comp.eps;']);
unix(['cd ~/Temp/painted_layer/; epsflatten comp.eps;']);
unix(['cd ~/Temp/painted_layer/; cp comp_flat.' destination(end-2:end) ' ' destination]);

