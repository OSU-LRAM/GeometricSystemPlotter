%Get box names and values from the working column    
function [box_names, box_active, box_values, box_enabled, plot_types,...
	plot_subtypes,merged_plot_subtypes, plot_style,xyplot]...
	= get_box_values(source_number_text,handles)

    %%%%
    %build the list of plot boxes for the column
    
    %classes of plot
    plot_types = {'vfield','CCF','bvi','disp','beta','dbeta','xy','xyopt'};
    
    %plot or subplot for each coordinate
    plot_style = {'plot','plot','subplot','subplot','plot','plot','single','single'};
    
    %text for the different components
    plot_subtypes = cat(1,repmat({{'X','Y','T'}},[6 1]),repmat({{'traj','net','BVI','cBVI'}},[2 1]));
	
	%Raw or optimized versions
	connection_types = {'','opt'};
	
	%merge plot subtypes
	for k = 1:length(plot_subtypes)
		for i = 1:length(connection_types)
			for j = 1:length(plot_subtypes{k})
				merged_plot_subtypes{k}{(i-1)*length(plot_subtypes{k})+j}...
					= [plot_subtypes{k}{j} connection_types{i}];
			end
		end
	end
	
	% Trim out redundant names from the beta and dbeta boxes, which don't
	% have opt versions
	merged_plot_subtypes{strcmp('beta',plot_types)}...
		(cellfun(@(x) ~isempty(x),(strfind(...
		merged_plot_subtypes{strcmp('beta',plot_types)},'opt')))) = [];
	
	merged_plot_subtypes{strcmp('dbeta',plot_types)}...
		(cellfun(@(x) ~isempty(x),(strfind(...
		merged_plot_subtypes{strcmp('dbeta',plot_types)},'opt')))) = [];
	
	% Trim out redundant names from the xy plot boxes (which handle
	% optimized and not optimized more separately
	merged_plot_subtypes{strcmp('xy',plot_types)}...
		(cellfun(@(x) ~isempty(x),(strfind(...
		merged_plot_subtypes{strcmp('xy',plot_types)},'opt')))) = [];
	
	merged_plot_subtypes{strcmp('xyopt',plot_types)}...
		(cellfun(@(x) ~isempty(x),(strfind(...
		merged_plot_subtypes{strcmp('xyopt',plot_types)},'opt')))) = [];

	
	%prime the arrays of box names and values
    box_names = cell(1,length(plot_types));%(2*length(plot_subtypes)+1)
    box_values = cell(size(box_names));
    box_enabled = cell(size(box_names));
    
    %build the list
    for i = 1:length(plot_types)
		
		box_names{i} = cell(1+length(merged_plot_subtypes{i}),1);
		box_values{i} = zeros(size(box_names{i}));
		box_enabled{i} = zeros(size(box_names{i}));                

		%category name
        box_names{i}{1} = plot_types{i};
        
        %category value
        box_values{i}(1) = ...
                get(handles.([plot_types{i} 'checkbox' source_number_text]),'Value');
            
        %category enabled state
        box_enabled{i}(1) = ...
				strcmp('on',get(handles.([plot_types{i} 'checkbox' source_number_text]),'Enable'));
            

        %specific plot names
        for j = 1:length(merged_plot_subtypes{i})
			
			%name
			box_names{i}{1+j} = ...
				[plot_types{i} plot_subtypes{j}];

			%value
			box_values{i}(1+j) = ...
				get(handles.([plot_types{i} merged_plot_subtypes{i}{j} 'checkbox' source_number_text]),'Value');

			%enabled state
			box_enabled{i}(1+j) = ...
				strcmp('on',get(handles.([plot_types{i} merged_plot_subtypes{i}{j} 'checkbox' source_number_text]),'Enable'));
			
        end
        
    end
    
    %Active boxes are both checked and enabled
    box_active = cellfun(@(x,y) x & y,box_values , box_enabled,'UniformOutput',false);
	

    
end