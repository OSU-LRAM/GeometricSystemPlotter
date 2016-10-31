function plots_to_make = initialize_plot_windows(box_active,plot_types,merged_plot_subtypes...
	,plot_style,hfuntype,stretchstate,stretchpath,handles,source_number_text)

	%%%%%
    %Determine how many windows to create
    
    plot_count = 0;                         %how many plots
    subplot_count = [];                     %how many subplots
    plots_to_make = struct([]);    %structure to hold all plot information
    %Loop through the categories
    for i = 1:length(plot_types)
        
        %only make plots for a category if that category is enabled and has
        %plots selected
        if box_active{i}(1) && any(box_active{i}(2:end))
			            
            %set the plots_to_make structure category for this plot set
            plots_to_make(end+1,1).category = plot_types{i}; %#ok<AGROW>
            
			%Appearance changing logic
			switch get(handles.(['contour' source_number_text]),'Value')

				case 0

					plots_to_make(end,1).style = 'surface';

				case 1

					plots_to_make(end,1).style = 'contour';

				otherwise

					error('Unexpected value for height function appearance toggle')

			end
			
			switch get(handles.(['color' source_number_text]),'Value')

				case 0
					
					plots_to_make(end,1).linestyle = 'mono';
					
				case 1
					
					plots_to_make(end,1).linestyle = 'cycle';
					
				otherwise
				
					error('Unexpected value for line color appearance toggle');
			end
				
			%Decide which kind of height function to show
			plots_to_make(end,1).hfuntype = hfuntype;

			% set the stretch on each component
			plots_to_make(end,1).stretch = ~strcmp(stretchpath,'null');
			plots_to_make(end,1).stretchpath = stretchpath;

			%set the plots_to_make components
            plots_to_make(end,1).components = merged_plot_subtypes{i}(box_active{i}(2:end));
			
			% Wrap an extra cell layer around the xy plots, to appear as
			% "one component"
			if strncmp(plots_to_make(end,1).category,'xy',2)
				plots_to_make(end,1).components = {plots_to_make(end,1).components};
			end
            
            %prime plot and subplot fields
            plots_to_make(end,1).plot = [];
            plots_to_make(end,1).subplot = [];
            
            %decide whether to make plots or subplots for the category
            switch plot_style{i}
                
                %if there's a plot for each component, add that many plots,
                %and note that there's one subplot for each of them
                case 'plot'
                    
                    %set the plot that will be associated with each
                    %component
                    for j = plot_count + (1:sum(box_active{i}(2:end)));
                        
                        plots_to_make(end,1).plot = [plots_to_make(end,1).plot;j];
                        plots_to_make(end,1).subplot = [plots_to_make(end,1).subplot; 1];
                        
                    end
                
                    plot_count = plot_count + sum(box_active{i}(2:end));
                    subplot_count = [subplot_count;ones(sum(box_active{i}(2:end)),1)]; %#ok<AGROW>
                    
                %if there's a subplot for each component, add one plot, and
                %note how many subplots
                case 'subplot'
                    
					for j = 1:2
					
						%set the plot that will be associated with all
						%components
						for k = 1:sum(box_active{i}((2:4)+(j-1)*3));

							plots_to_make(end,1).plot = [plots_to_make(end,1).plot; plot_count+1];
							plots_to_make(end,1).subplot = [plots_to_make(end,1).subplot; k];

						end

						if (sum(box_active{i}((2:4)+(j-1)*3)) > 0)
							plot_count = plot_count + 1;
							subplot_count = [subplot_count;sum(box_active{i}((2:4)+(j-1)*3))]; %#ok<AGROW>
						end
                    
					end
					
				case 'single'
					
					% for the xy plots, which have a single window
					plots_to_make(end,1).plot = [plots_to_make(end,1).plot;plot_count+1];
					plots_to_make(end,1).subplot = [plots_to_make(end,1).subplot; 1];
					
					plot_count = plot_count + 1;
					subplot_count = [subplot_count;1];
					
                otherwise
                    
                    error('plot style is not ''plot'' nor ''subplot'' nor ''single'' ')
                    
            end
            
        end
        
    end
    
    
    %%%%%%
    %Generate the plot axes
    
    %make more columns than rows
    n_columns = ceil(sqrt(plot_count));
    n_rows = ceil(plot_count/n_columns);
    
    %prime a cell array of axes handles to plot into
    thumbnail_axes = cell(plot_count,1);
    
    
    %generate the plot axes
    for i = 1:plot_count
        
        %Prime the array of thumbnail axes for this plot
        thumbnail_axes{i} = zeros(1,subplot_count(i));
        
        %get the row/column position of the current plot - note that
        %subplot indexing goes across rows first, so we need [J,I] instead
        %of [I,J]
        [J,I] = ind2sub([n_columns,n_rows],i);
        
        for j = 1:subplot_count(i)
            
            %convert to row/column position with finer subgroup
            I_f = subplot_count(i)*(I-1)+j;
            J_f = J;
            
            %Get index from fine row/column position, remembering to use
            %the transpose
            IND = sub2ind([n_columns,n_rows*subplot_count(i)],J_f,I_f);
            
            thumbnail_axes{i}(j) = ...
                subplot(n_rows*subplot_count(i),n_columns,IND,'Parent',handles.plot_thumbnails);
            
        end
        
    end
            
        
    %%%%%%%%%%
    %Associate axes handles with plot commands
    
    
    %Loop over all plot categories that will be used
	for i = 1:length(plots_to_make)
        
        %loop over all components
        for j = 1:length(plots_to_make(i).components)
            
            %get the axes to plot to
            plots_to_make(i).axes(j,1) = thumbnail_axes{plots_to_make(i).plot(j)}(plots_to_make(i).subplot(j));
            
        end
        
	end
	
end