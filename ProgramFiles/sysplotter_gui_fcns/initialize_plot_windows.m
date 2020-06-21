% function plots_to_make = initialize_plot_windows(box_active,plot_types,merged_plot_subtypes...
% 	,plot_style,CCFtype,stretchstate,handles,source_number_text)

function plots_to_make = initialize_plot_windows(box_names,...
                                                 box_active,...
                                                 plot_style,...
                                                 stretchstate,...
                                                 current_stretch,...
                                                 plot_types,...
                                                 plot_subtypes,...
                                                 plot_coordinates,...
                                                 handles,...
                                                 source_number_text,...
                                                 CCFtype)

	%%%%%
    %Determine how many windows to create
    
    plot_count = 0;                         %how many plots
    subplot_count = [];                     %how many subplots
    plots_to_make = struct('components',{},...    %structure to hold all plot information
                           'plot',[],...
                           'subplot',[],...
                           'category',[],...
                           'linestyle',[],...
                           'style',[],...
                           'CCFtype',[],...
                           'stretch',[],...
                           'axes',[],...
                           'stretch_name',[]);
    
    
    
    % Loop through the checkboxes
    
    % Iterate through different kinds of plots
    for idx = 1:numel(box_names)
        
        
        % Iterate through coordinate choices
        for idx2 = 1:numel(box_names{idx})
    
                % Check if there are any plots of this type
                plot_type_active = [box_active{idx}{idx2}{:}];

                

                
                if any(plot_type_active)


                    %decide whether to make plots or subplots for the category
                    switch plot_style{idx}

                        %if there's a plot for each component, add that many plots,
                        %and note that there's one subplot for each of them
                        case 'plot'



                                %set the plot that will be associated with each
                                %component

                                % iterate over the plots in the category
                                for idx3 = 1:numel(box_active{idx}{idx2}) 

                                    plot_subtype_active = box_active{idx}{idx2}{idx3};

                                    if plot_subtype_active
                                        
                                        % increment the plot counter
                                        plot_count = plot_count+1;

                                        %fill in plot and subplot fields
                                        plots_to_make(end+1,1).plot = plot_count;
                                        plots_to_make(end,1).subplot = 1;
                                        %set the plots_to_make components

                                        % Name this component
                                        plots_to_make(end,1).components =...
                                            {[plot_subtypes{idx}{idx3} plot_coordinates{idx}{idx2}]};

                                        % Note that there is one subplot in
                                        % this plot
                                        subplot_count = [subplot_count;1];
                                        
                                        plots_to_make(end,1) = add_category_info(plots_to_make(end,1),...
                                             plot_types,...
                                             idx,...
                                             handles,...
                                             CCFtype,...
                                             stretchstate,...
                                             current_stretch,...
                                             source_number_text);

    % 
    %                                      % Iterate through subplots
    %                                      for idx3 = 1:numel(box_names{idx}{idx2})
    % 
    %                                         plot_subtype_active = box_active{idx}{idx2}{idx3};
    % 
    %                                         if plot_subtype_active
    % 
    %                                             plots_to_make(end,1).components = ...
    %                                                 [plots_to_make(end,1).components,... 
    %                                                  [plot_subtypes{idx}{idx3} plot_coordinates{idx}{idx2}]...
    %                                                  ]
    % 
    %                                         end
    % 
    %                                     end
    % 
    %                                     % increment the plot counter
    %                                     plot_count = plot_count+1;
    % 
    %                                     % Add a 1-subplot entry to the subplot counter
    %                                     subplot_count = [subplot_count; 1];
    % 
    %                                     % Label the plot and subplot in the structure
    %                                     plots_to_make(end,1).plot = plot_count;
    %                                     plots_to_make(end,1).subplot = [plots_to_make(end,1).subplot; 1];

                                
                                    end
                                    
                                end
                                
                        %if there's a subplot for each component, add one plot, and
                        %note how many subplots
                        case 'subplot'

                            % increment the plot counter
                            plot_count = plot_count+1;

                            % iterate over the plots in the category
                            %for idx3 = numel(plots_to_make(idx).components)
                            for idx3 = 1:numel(box_active{idx}{idx2}) 


                                plot_subtype_active = box_active{idx}{idx2}{idx3};

                                if plot_subtype_active
                                    
                                    % Label the plot and subplot in the structure
                                    plots_to_make(end+1,1).plot = plot_count;
                                    plots_to_make(end,1).subplot = [plots_to_make(end,1).subplot; idx3];


                                    plots_to_make(end,1).components = ...
                                        [plots_to_make(end,1).components,... 
                                         {[plot_subtypes{idx}{idx3} plot_coordinates{idx}{idx2}]}...
                                         ];
                                     
                                    plots_to_make(end,1) = add_category_info(plots_to_make(end,1),...
                                     plot_types,...
                                     idx,...
                                     handles,...
                                     CCFtype,...
                                     stretchstate,...
                                     current_stretch,...
                                     source_number_text);

                                end



                             end

                            % Tell the subplot counter how many subplots appear
                            % in this plot
                            subplot_count = [subplot_count; idx3];

                            
                            
                        % for plots which have a single item    
                        case 'single'

                            % increment the plot counter
                            plot_count = plot_count+1;
                            
                            plots_to_make(end+1,1).plot = plot_count;
                            plots_to_make(end,1).subplot = [plots_to_make(end,1).subplot; 1];
                            
                            % iterate over the plots in the category
                            %for idx3 = numel(plots_to_make(idx).components)
                            for idx3 = 1:numel(box_active{idx}{idx2}) 
                                
                                
                                    plots_to_make(end,1) = add_category_info(plots_to_make(end,1),...
                                     plot_types,...
                                     idx,...
                                     handles,...
                                     CCFtype,...
                                     stretchstate,...
                                     current_stretch,...
                                     source_number_text);

                                plot_subtype_active = box_active{idx}{idx2}{idx3};

                                if plot_subtype_active
                                    
                                    % Label the plot and subplot in the structure
                                    
                               


                                    plots_to_make(end,1).components = ...
                                        [plots_to_make(end,1).components,... 
                                         {[plot_subtypes{idx}{idx3} plot_coordinates{idx}{idx2}]}...
                                         ];
                                     


                                end



                            end
                             
                            plots_to_make(end,1).components = {plots_to_make(end,1).components};

                            % Tell the subplot counter how many subplots appear
                            % in this plot
                            subplot_count = [subplot_count; 1];                            
                            
%                             plots_to_make(end+1,1).plot = plot_count;
%                             plots_to_make(end,1).subplot = 1;
% 
%                             plots_to_make(end,1).components = ...
%                                 [plots_to_make(end,1).components,... 
%                                  {[plot_subtypes{idx}{1} plot_coordinates{idx}{idx2}]}...
%                                  ];
%                             
%                             subplot_count = [subplot_count;1];
% 
%                             plots_to_make(end,1) = add_category_info(plots_to_make(end,1),...
%                                      plot_types,...
%                                      idx,...
%                                      handles,...
%                                      CCFtype,...
%                                      stretchstate,...
%                                      source_number_text);                            

                    
                        otherwise

                            error('plot style is not ''plot'' nor ''subplot'' nor ''single'' ')

                end
                        
%                         % If there are any active boxes of this type, add them to the
%             % plots_to_make structure 
% 
%             plot_type_active = [box_active{idx}{idx2}{:}];
% 
%             if any(plot_type_active)



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
        
        plots_to_make(i).axes = zeros(length(plots_to_make(i).subplot),1);
        
        %loop over all components
        for j = 1:length(plots_to_make(i).subplot)
            
            %get the axes to plot to
            plots_to_make(i).axes(j,1) = thumbnail_axes{plots_to_make(i).plot(j)}(plots_to_make(i).subplot(j));
            
        end
        
	end
	
end


function plots_to_make_i = add_category_info(plots_to_make_i,...
                                             plot_types,...
                                             idx,...
                                             handles,...
                                             CCFtype,...
                                             stretchstate,...
                                             current_stretch,...
                                             source_number_text)

    %set the plots_to_make structure category for this plot type
    plots_to_make_i.category = plot_types{idx};

    % Decide if scalar functions of two variables should be plotted
    % as contours or surfaces 
    switch get(handles.(['contour' source_number_text]),'Value')

        case 0

            plots_to_make_i.style = 'surface';

        case 1

            plots_to_make_i.style = 'contour';

        otherwise

            error('Unexpected value for CCF appearance toggle')

    end
    
%     if strcmp(current_stretch,'metric_surface')
%                 
%         plots_to_make_i.style = 'pcolor';
%                 
%     end
    % Decide if monocolor or a cycling color set should be used for
    % multi-line plots
    switch get(handles.(['color' source_number_text]),'Value')

        case 0

            plots_to_make_i.linestyle = 'mono';

        case 1

            plots_to_make_i.linestyle = 'cycle';

        otherwise

            error('Unexpected value for line color appearance toggle');
    end

    %Decide which kind of constraint curvature function to show
    plots_to_make_i.CCFtype = CCFtype;

    % set the stretch on each component. Menu item 1 is no stretch,
    % 2 is stretch (and if we offer multiple stretches in the
    % future, this will pass on higher values into the plotting
    % function)
    plots_to_make_i.stretch = stretchstate-1;
    plots_to_make_i.stretch_name = current_stretch;
                
end