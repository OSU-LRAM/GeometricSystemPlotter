function painter_layer2(f,layers,imageformat,destination,printoptions)
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
	
	
	
% 	% Set the renderer
% 	%check the elements of the layer for surfaces
% 	surf_found = strcmp('surface',get(layers{i},'type'));
% 	if any(surf_found)
% 		renderer_string = {'-painters'};
% 	else
% 		renderer_string = {'-painters'};
% 	end
		
	% Set the image print driver string
	switch imageformat{i}
		
		case 'vector'
			
			print_driver = '-depsc2';
			
		case 'raster'
			
			print_driver = '-dpng';
			
			% greenscreen the background
			oldbackgroundcolor = get(f,'Color');
			set(f,'Color',[0 1 0]);
			
		otherwise
			
			error('Unsupported image format specified')
			
	end
	
	%print this layer
	print(f,print_driver,'-loose','-painters',printoptions{:},['~/Temp/painted_layer/' num2str(i)])		
	%epsfontsubs(f,['~/Temp/painted_layer/' num2str(i)],printoptions{:})
	
	% Convert to pdf
	switch imageformat{i}
		
		case 'vector'
			
			% Get the bounding box, flatten the file, and then restore the
			% bounding box
			
			
			%%%%%%%%%%%
			% get bounding boxes
			dirname = '~/Temp/painted_layer/';
			filename = [num2str(i) '.eps'];
			fid = fopen(fullfile(dirname,filename));

			% skip the first string
			% loop down strings in file
			j = 0;
			while ~feof(fid)

				%keep count of where we are
				j = j+1;

				% read the string
				s = fgetl(fid);

				if strncmp(s,'%%Bounding',10)
					bboxstr = strrep(s,'%%BoundingBox: ','');

					% store bounding box info into array
					bbox = str2num(bboxstr);

				end

				if strncmp(s,'%%HiResBounding',15)
					hiresbboxstr = strrep(s,'%%HiResBoundingBox: ','');

					% store bounding box info into array
					hiresbbox = str2num(hiresbboxstr);

				end

			end

			% close the file
			fclose(fid);
			
			%%%%%%%%%
			% Flatten the file
			unix(['cd ~/Temp/painted_layer/; epsflatten ' num2str(i) '.eps;']);
			
			%%%%%%%
			% Restore the bounding box
			fido = fopen(fullfile(dirname,'bboxtemp'),'w'); % output temp file

			fid = fopen(fullfile(dirname,filename));

			% loop down strings in file
			j = 0;
			while ~feof(fid)

				%keep count of where we are
				j = j+1;

				% read the string
				s = fgetl(fid);

				if strncmp(s,'%%Bounding',10)
					s = ['%%BoundingBox: ' num2str(bbox)];
				end

				if strncmp(s,'%%HiResBounding',15)
					s = ['%%HiResBoundingBox: ' num2str(hiresbbox)];
				end

				fprintf(fido,'%s\n',s);


			end

			% close the file
			fclose(fid);
			fclose(fido);

			system(['mv ', fullfile(dirname,'bboxtemp'), ' ',fullfile(dirname,filename)]);

			%%%%%%%%%%
			% move the flattened pdf into position for composition
			unix(['cd ~/Temp/painted_layer/; epstopdf ' num2str(i) '.eps; mv ' num2str(i) '_flat.pdf ' num2str(i) '.pdf']);
			
		case 'raster'
			
 			unix(['cd ~/Temp/painted_layer/; convert ' num2str(i) '.png -transparent \#0f0 ' num2str(i) '.pdf']);
%			unix(['cd ~/Temp/painted_layer/; greenscreentopdf ' num2str(i) '.png']);
			
			set(f,'Color',oldbackgroundcolor);
			
		otherwise
			
			error('Unsupported image format specified')
			
	end


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

% Start by putting layer 2 on layer 1
unix(['cd ~/Temp/painted_layer/; pdftk 2.pdf background 1.pdf output combinedfile.pdf'])

% loop over the rest of the layers
for i = 3:length(layers)
	unix(['cd ~/Temp/painted_layer/; mv combinedfile.pdf oldcombinedfile.pdf; pdftk ' num2str(i) '.pdf background oldcombinedfile.pdf output combinedfile.pdf']);
end

%%%%%%%
% send file to destination
switch destination(end-2:end)
	
	case 'pdf'
		
		unix(['cd ~/Temp/painted_layer/; cp combinedfile.pdf ' destination]);
		
	case 'eps'
		
		unix(['cd ~/Temp/painted_layer/; pdftoeps combinedfile.pdf; cp combinedfile.eps ' destination]);

	otherwise
		
		error('unsupported final export format')
		
end
		


% %build up the arguments for epscompose
% filenames = [num2str((1:length(layers))'),repmat('.eps 0 0 1 0 ',length(layers),1)]';
% 
% unix(['cd ~/Temp/painted_layer/; epscompose ' (filenames(:))' ' > comp.eps;']);
% unix(['cd ~/Temp/painted_layer/; epsflatten comp.eps;']);
% unix(['cd ~/Temp/painted_layer/; cp comp_flat.' destination(end-2:end) ' ' destination]);
% 
