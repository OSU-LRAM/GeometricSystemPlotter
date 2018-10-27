function illustrate_shapespace(system,target)

    geometry = system.geometry;
    visual = system.visual;

    load('sysplotter_config','Colorset')


    % Specify plot number as 177 unless specified otherwise
    if ~exist('target','var')
        target = 177;
    end

    % Take the display grid from the provided display structure. Make sure that
    % the components of the grid are in a column cell array to avoid any
    % row/column problems below
    paramgrid = visual.grid(:);
    n_dim=numel(visual.grid);
    test2=cell(1,n_dim-1);
    test3=cell(1,n_dim-2);
    test2{1}=1;
    for i=1:length(test3)
        test2{1,i+1}=1;
        test3{1,i}=1;
    end


    % Get the spacing of parameter values in the first two dimensions
    [~,grid_spacing_x] = gradient(paramgrid{1});
    [grid_spacing_y,~] = gradient(paramgrid{2});

    % Create a figure if necessary
    try 
        target_type = get(target,'type');
    catch
        target_type = 'DNE';
    end

    if strcmp(target_type,'axes')
        axh = target;
    else
        fh = figure(target);
        clf(fh);
        fh = figure(target);
        set(fh,'name',['Baseframe: ' geometry.baseframe]);

        axh = axes('Parent',fh);
    end

    % Set up axes
    axis(axh,'equal')
    hold(axh,'on')
    % Set tick spacing
    set(axh,'XTick', paramgrid{1}(:,test2{:}));
    set(axh,'YTick', paramgrid{2}(1,:,test3{:}));
    set(axh,'XLim', [min(paramgrid{1}(:))-grid_spacing_x(1)/2,max(paramgrid{1}(:))+grid_spacing_x(end)/2]);
    set(axh,'YLim', [min(paramgrid{2}(:))-grid_spacing_y(1)/2,max(paramgrid{2}(:))+grid_spacing_y(end)/2]);
    % Figure stylistic settings
    set(axh, 'Box','on'); % Repeat tickmarks and axis lines on the top and left sides
    grid(axh,'on'); % show a grid at the specified shape points

    % The system length scale should be as long as the minimum spacing between nearby
    % elements
    blnth = .9 * min([grid_spacing_x(:);grid_spacing_y(:)]);

    % Scale the geometry defined in the system file
     geometry.length = blnth;


    % Specfication of face colors for different parts of the system. Currently
    % support primary portion of the system in white with a black outline, and
    % secondary portion in the spot color with a black outline, and a third
    % portion in the secondary color
    fillcolors = {'w',Colorset.spot,Colorset.secondary};

    % define a color list
    colorlist = {Colorset.spot, Colorset.secondary, 'k'};
    
    % In 4-dimensional systems, we need to select colors based on the value
    % of the 4th shape variable.  
    if n_dim==4
        maxp=max(max(max(max(paramgrid{1}))));
        minp=min(min(min(min(paramgrid{1}))));
    end


    for idx = 1:numel(paramgrid{1})

        % Extract the parameter values for this combination of parameters
        p = zeros(size(paramgrid));
        for idx2 = 1:numel(paramgrid)
            p(idx2) = paramgrid{idx2}(idx);
        end
        
        % Code to make sure that for systems with 4 shape variables, the
        % swimmer with 4th shape variable zero is plotted in black and
        % plotted last for a given spot on alpha1-alpha2-alpha3 space.
        if n_dim==4
            value1=(length(paramgrid{4}(:,1)))^3;
            p(4)=paramgrid{4}(mod(idx+2*value1-1,numel(paramgrid{1}))+1);
        end

        %Generate the backbone locus
        B = generate_locomotor_locus(geometry,p,visual);

        % Force B to be a column cell array
        if ~iscell(B)
            B = {B};
        end
        B = B(:);

        for idx2 = 1:numel(B)

            % Iterate over each piece within the body

            % Force B{idx2} to be a column cell array
            if ~iscell(B{idx2})
                B{idx2} = {B{idx2}};
            end
            B{idx2} = B{idx2}(:);

            for idx3 = 1:numel(B{idx2})
                % Taking the coordinates of the swimmer in R^2 and embedding it in the shape space 
                if n_dim==3
                    B{idx2}{idx3}(3,:)=zeros(1,length(B{idx2}{idx3}(1,:)));
                end
                if n_dim==4
                    B{idx2}{idx3}(3,:)=zeros(1,length(B{idx2}{idx3}(1,:)));
                    B{idx2}{idx3}(4,:)=zeros(1,length(B{idx2}{idx3}(1,:)));
                end
                % transform the body elements to their position on the grid
                if n_dim<4
                    for idx4 = 1:n_dim
                        B{idx2}{idx3}(idx4,:) = B{idx2}{idx3}(idx4,:) + p(idx4);
                    end
                else
                    for idx4 = 1:3
                        B{idx2}{idx3}(idx4,:) = B{idx2}{idx3}(idx4,:) + p(idx4);
                    end
                end
                % draw the backbone at the specified location
                if n_dim==2
                    bp{idx,idx2,idx3} = patch('Xdata',B{idx2}{idx3}(1,:),'Ydata',B{idx2}{idx3}(2,:),'Parent',axh,'EdgeColor','k','FaceColor',fillcolors{idx2});
                end
                if n_dim==3
                    bp{idx,idx2,idx3} = patch('Xdata',B{idx2}{idx3}(1,:),'Ydata',B{idx2}{idx3}(2,:),'Zdata',B{idx2}{idx3}(3,:),'Parent',axh,'EdgeColor','k','FaceColor',fillcolors{idx2});
                end
                if n_dim==4
                    if p(4)==0
                        bp{idx,idx2,idx3} = patch('Xdata',B{idx2}{idx3}(1,:),'Ydata',B{idx2}{idx3}(2,:),'Zdata',B{idx2}{idx3}(3,:),'Parent',axh,'EdgeColor','k');
                    elseif p(4)==minp
                        bp{idx,idx2,idx3} = patch('Xdata',B{idx2}{idx3}(1,:),'Ydata',B{idx2}{idx3}(2,:),'Zdata',B{idx2}{idx3}(3,:),'Parent',axh,'EdgeColor',Colorset.spot);
                    else 
                        bp{idx,idx2,idx3} = patch('Xdata',B{idx2}{idx3}(1,:),'Ydata',B{idx2}{idx3}(2,:),'Zdata',B{idx2}{idx3}(3,:),'Parent',axh,'EdgeColor',Colorset.secondary);
                    end
                end
            end
        end
    end
end