function [convert,sol] = fast_flatten_metric(grid,metric,mask)
% Calculate a change of coordinates to optimally flatten a metric tensor.
% This function uses a linear-distortion metric (related to Tissot's
% indicatrix) and finds the optimal parameterization through a
% spring-relaxation method
%
% Inputs:
% grid is a cell array of the output of ndgrid
% metric is a 2x2 cell array that contains the tensor components evaluated on the grid
% mask is an optional function that evaluates to true for all node locations inside a
%       region of interest and false elsewhere
%         (e.g., "@(x,y) x^2+y^2 <= 1" removes the influence of any grid
%         points outside the unit circle)
%       
%
% Outputs:
% convert: a structure of function handles:
%	old_to_new_points: returns (x_new,y_new) corresponding to input (x_old,y_old)
%	new_to_old_points: returns (x_old,y_old) corresponding to input (x_new,y_new)
%	jacobian: returns the 2x2 jacobian from old to new coordinates at (x_old,y_old)
%	new_metric: returns the 2x2 metric in the new coordinates, at (x_old,y_old)
% sol: the ODE solution describing the deformation path of the nodal points
%  in the parameterization from their start to final values
% sol_interp: Converts old to new points (in the manner of
% convert.old_to_new_points) that also takes in a time t (from 0 to 1)
% during the relaxation period, to say how far along the evolution this
% point conversion should be made

	
    %Build springs from grid
	[springs, blocks] = generate_springs(grid{1});

	% Find the initial lengths of the springs (the initial grid separations
	[start_lengths,~,start_deltas] = get_spring_lengths_and_azimuths(springs,grid{:});

	% Calculate the lengths the springs would be if the metric could be
	% flattened entirely
	[neutral_lengths,mean_neutral_length] = get_spring_neutral_lengths(springs,blocks,start_deltas,grid{:},metric);
    
    % Scale the initial positions by the ratio between the mean neutral
    % length and the mean initial length
    x_scaled = (grid{1}/geomean(start_lengths))*mean_neutral_length;
    y_scaled = (grid{2}/geomean(start_lengths))*mean_neutral_length;

    
    %%%%%
    % Masking functions for non-rectangular regions of the shape space
    
    % If no masking function was specified, make the trivial mask
    if ~exist('mask','var')
        mask = @(x,y) ones(size(x));
    end
    
    % Apply the mask to the points
    masked_points = mask(grid{1}(:),grid{2}(:));
    
    % Generate a mask for the springs, that returns zero if either end of
    % the spring is outside the target region.
    masked_springs = min(masked_points(springs(:,1:2)),[],2);
    
    % Append the masked springs attribute onto the end of the spring array
    springs = [springs,masked_springs]; 
    
    
    %%%%%%%%%%
    % Processing step
    %Relax the springs
	[final_x,final_y,sol] = relax_springs(x_scaled,y_scaled,springs,neutral_lengths,0.01);
    
	% Convert the metric to each location's tensor at a single
	% location
	metric = celltensorconvert(metric);
	
		
	
	%%%%%%%%%%%%
	%%%%%%%%%%%%
	%Build the output functions
	
	% Convert points
	convert.old_to_new_points = @(x_old,y_old) convert_points(grid{:},final_x,final_y,x_old,y_old);
	Fx = TriScatteredInterp([final_x(:) final_y(:)],grid{1}(:));
	Fy = TriScatteredInterp([final_x(:) final_y(:)],grid{2}(:));
	convert.new_to_old_points = @(x_new,y_new) multiTriInterp(Fx,Fy,x_new,y_new);
	
	% Jacobian from old to new tangent vectors
	final_jacobian = find_jacobian(grid{:},final_x,final_y);
	convert.jacobian = @(x_p,y_p) interpolate_cellwise_tensor(grid{:},x_p,y_p,final_jacobian);
	
	% Metric in new space
	final_new_metric = cellfun(@(j,m) j'\m/j,celltensorconvert(final_jacobian),(metric),'UniformOutput',false);
	convert.new_metric = @(x_p,y_p) interpolate_cellwise_tensor(grid{:},x_p,y_p,celltensorconvert(final_new_metric));
	convert.old_metric = @(x_p,y_p) interpolate_cellwise_tensor(grid{:},x_p,y_p,celltensorconvert(metric));
    
    
	
end