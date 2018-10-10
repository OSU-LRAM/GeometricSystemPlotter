function [gradE, E]...
	= helmholtz(grid,V,weight,data_storage_location,grid_regularity)

% Discrete Hodge-Helmholtz-esqe decomposition on a rectangular region, with
% optimization to store the transform matrix. 
%
% Inputs are:
%  grid: cell of ndgridoutputs
%  V: vector field structure -- cell array of arrays, nth cell
%    corresponding to nth component
%  weight: relative importance of each point (optional, defaults to even weighting)
%  data_storage_location: directory in which to cache results of
%    calculation (optional, defaults to the current directory)
%  grid_regularity: future option, for specifying non-regular grids.
%    (optional, defaults to only valid option, 'even grid')
%
% Outputs are: 
%  gradE: Curl-free component of the vector field
%  E: Scalar potential function whose gradient is the curl-free component,
%  and which has a value of zero at the top left corner
%

	% Set default arguments
	if ~exist('weight','var')
		weight = ones(size(grid{1}));
	end
	
	if ~exist('data_storage_location','var')
		data_storage_location = '.';
	end

	if ~exist('grid_regularity','var')
		grid_regularity = 'even grid';
	end
	
	%Get the matrix mapping from velocities to projected scalar potential	
	[Sinv,nodes,cubes] = preprocessor(grid,weight,data_storage_location,grid_regularity); %#ok<NASGU,ASGLU>
	
	%Multiply Ainv by the input vector field to get the scalar potential
	%function solution that has E_1 = 0 (as specified in the preprocessor)
	[E] = processor(V,Sinv);


	%Process E_test to get the desired outputs
	[gradE] = postprocessor(grid,E);
	
end
	

function [Sinv,nodes,cubes] = preprocessor(grid,weight,data_storage_location,grid_regularity)
	
	% get the location of this file
	[pathstr] = fileparts(which([mfilename '.m']));

	%%%%%%%%%
	%First, check if we've already calculated this matrix

	recalc = 1; %default value of flag for recalculation
		
	%Name the directory in which we store precalculated data
	datadir = fullfile(data_storage_location, [mfilename '_data']);
	
	%If the precalc data directory doesn't exist, make it
	if ~exist(datadir,'dir')
		
		mkdir(datadir);
		
	end
	
	%Hash the input values, and make a string from the return
	gridcolumns = grid_to_columns(grid);
	input_params = [range(gridcolumns,1) weight(:)'];
	input_hash = num2str(hash(input_params,'MD5'));

	%Compare the date of the hashnamed data file with the date of this
	%function (getfiledate will do the right thing if the file doesn't
	%exist yet).
	D1 = dir(pathstr); %list of files in this directory
	D2 = dir(datadir); %list of precalculated files
	
	%get the file dates
	helmholtz_date = getfiledate(D1,[mfilename '.m']);
	data_date = getfiledate(D2,[input_hash '.mat']);
	
	%If the data file has been updated more recently than the
	%helmholtz function, don't need to rerun the calculation
	if (data_date > helmholtz_date)
		
		%Also, check the actual size and spacing to protect against hash
		%collisions
		load([datadir '/' input_hash '.mat'],'grid_params')
		params_test = sum(grid_params == input_params); %#ok<NODEF>
		if params_test == length(grid_params)
			
			recalc = 0;
			
		end
		
	end

	
	%%%%%%%%%%%%%
	%Decide whether or not to recalculate
	
	if recalc
	
	
		%%%%%%%%
		%Construct a cuboid mesh
		
		% Get the corresponding normalized mesh
		[nodes,cubes] = hypercube_mesh(grid);

	
		
		%Separate out useful values from grid size
		N = size(gridcolumns,1); %number of nodes
		n_dim = size(gridcolumns,2); % number of dimensions over which fields are defined

		%Prime the S1 and S2 matrices
		S1 = zeros(N);
		S2 = zeros(N,n_dim*N);
		
		% Define the (n-)linear shape functions for the finite elements
		[shape_functions, shape_dfunctions] = hypercube_element_shape_functions(n_dim);
		
		% Extract the number of shape functions
		n_shape = numel(shape_functions);
		
		% Group the shape functions
		shape_functions_group = @(X) cellfun(@(f) f(X), shape_functions,'UniformOutput',false);
		
		% Group the shape function derivatives by row (each row contains
		% the derivatives of a given shape function with respect to all
		% dimensions)
		shape_dfunctions_group = cell(n_shape,1);
		for i = 1:n_shape
			shape_dfunctions_group{i,1} = @(X) cellfun(@(f) f(X), shape_dfunctions(i,:),'UniformOutput',false);
		end
		
		% Group the shape function derivatives by column (each column contains
		% the derivatives of a given shape function with respect to all
		% dimensions)
		shape_dfunctions_group_by_col = cell(1,n_dim);
		for i = 1:n_dim
			shape_dfunctions_group_by_col{1,i} = @(X) cellfun(@(f) f(X), shape_dfunctions(:,i),'UniformOutput',false);
		end
		
		% Get the quadrature points and weights for the elements
		[quad_points, quad_weights] = hypercube_quadrature_init(n_dim);
		
		% Calculate the shape functions and their derivatives at the
		% quadrature points (still in coordinates local to the elements)
		shape_at_quad = shape_functions_group(quad_points);
		dshape_at_quad = cellfun(@(f) f(quad_points),shape_dfunctions_group,'UniformOutput',false);

		%%%%%%%%%
		% Create the Jacobian from the -1 to 1 hypercube to a
		% given element

		% Interpolation functions (sum of the basis functions, weighted by
		% the locations of the actual node elements)
		Te = cell(n_dim,1);
		for k = 1:n_dim
			Te{k} = @(Eta,element) cell_weighted_sum(gridcolumns(cubes(element,:),k),shape_functions_group(Eta));
		end

		% Derivative interpolation functions (sum of the basis derivative
		% functions, weighted by the locations of the actual elements
		dTe = cell(n_dim,n_dim);
		for k = 1:size(dTe,1)
			for m = 1:size(dTe,2)
				dTe{k,m} = @(Eta,element) cell_weighted_sum(gridcolumns(cubes(element,:),k),shape_dfunctions_group_by_col{m}(Eta));
			end
		end

		% Build up the Jacobian function from dTe -- note that this expects
		% a single point input.
		J = @(Eta,element) cellfun(@(f) f(Eta,element),dTe);
				
		%%%%%
		%Build up the S1 and S2 matrices
		
		% Select grid regularity (
		switch grid_regularity
				
			case 'even grid'
				
				% Jacobian is the same for all grid elements and constant
				% over the cell
				J_element = J(num2cell(zeros(n_dim,1)),1);
				detJ_element = det(J_element);
				
				% basis function derivative in local coordinates for each
				% basis function
				gradPhi_local_quad = arrayfun(@(x) [dshape_at_quad{x}{:}]',1:size(cubes,2),'UniformOutput',false);

				% transform derivative into real coordinates
				gradPhi_transformed_quad = cell(size(gradPhi_local_quad));
				for i = 1:numel(gradPhi_local_quad)
					gradPhi_transformed_quad{i} = cell2mat(...
						cellfun(@(x) J_element\x,...      
						mat2cell(gradPhi_local_quad{i},size(gradPhi_local_quad{i},1)...
							,ones(1,size(gradPhi_local_quad{i},2)) ),'UniformOutput',false) );
				end

				%%%%%				
				% Fill in the dependency of the potential-function gradient
				% on the potential function at the nodes

				% Values to multiply by the potential function values
				% to get its gradient
				P_local_quad = arrayfun(@(c) [dshape_at_quad{c}{:}]',1:n_shape,'UniformOutput',false);


				% Transform the potential gradient into real
				% coordinates
				P_transformed_quad = cell(size(P_local_quad));
				for m = 1:numel(P_transformed_quad)

					P_transformed_quad{m} = cell2mat(...
						cellfun(@(x) J_element\x,...
						mat2cell(P_local_quad{m},size(P_local_quad{m},1)...
							,ones(1,size(P_local_quad{m},2)) ),'UniformOutput',false) );

				end
						
				% The vector fields are independent of the spacing, so they
				% can be interpolated using simply the shape function
				% values at the quadrature points
				R_quad = cellfun(@(x)x',shape_at_quad','UniformOutput',false);
						
				%%%%%%%%%%%%%%%%%%%%%%%
				% Get the dot product of the basis gradient with the
				% potential gradient and the input vector field

				% full dot product for the potential function
				gradPhi_dot_P_quad = cell(size(gradPhi_local_quad));
				for i = 1:numel(gradPhi_local_quad)
					gradPhi_dot_P_quad{i} ...
						= cellfun(@(p) sum(gradPhi_transformed_quad{i}.*p,1)...
							,P_transformed_quad,'UniformOutput',false);
				end

				% keep components separate for the input vector field
				gradPhi_dot_R_quad = cell(size(gradPhi_local_quad));
				for i = 1:numel(gradPhi_local_quad) 
					gradPhi_dot_R_quad{i} ...
						= cellfun(@(r) gradPhi_transformed_quad{i}.*repmat(r,[n_dim,1])...
							,R_quad,'UniformOutput',false);
				end

				% Sum the dot products over the quadrature weights to get
				% the integral value's relationship to the node values in
				% this cell, also multiplying by the cell's jacobian
				% determinant
				gradPhi_dot_P = cell(size(gradPhi_local_quad));
				gradPhi_dot_R = cell(size(gradPhi_local_quad));
				for i = 1:numel(gradPhi_local_quad) 

					gradPhi_dot_P{i} = cellfun(@(x) detJ_element*x*quad_weights(:),gradPhi_dot_P_quad{i});
					gradPhi_dot_R{i} = cell2mat(cellfun(@(x) detJ_element*x*quad_weights(:),gradPhi_dot_R_quad{i},'UniformOutput',false));
				
				end
				
				% Loop over each node, filling in the integral data from
				% the appropriate precomputed value
				for i = 1:N

					% Get the list of cubes containing this node
					cubes_at_node = find(any(cubes==i,2));
					cube_list = cubes(cubes_at_node,:); %#ok<FNDSB>
				
					% Loop over each cube attached to the node
					for j = 1:size(cube_list,1)

						% Get which corner (and thus basis function) of the cube is
						% the current node
						cornerI = cube_list(j,:)==i;

						% Slot the dotproduct values into a vector that pairs them
						% with the nodal values they correspond to
						S1(i,cube_list(j,:)) = S1(i,cube_list(j,:)) + gradPhi_dot_P{cornerI};

						% Offset each component's-worth of R values by the number
						% of nodes
						for m = 1:n_dim
							S2(i,cube_list(j,:)+(m-1)*N) = S2(i,cube_list(j,:)+(m-1)*N) + gradPhi_dot_R{cornerI}(m,:);
						end


					end
					
				end
				
			otherwise
				
				error('Unknown matrix construction algorithm')
				
		end
				
		% Multiply the matrices by the node weighting factors
		%S1 = diag(weight(:))*S1;
		S2 = diag(weight(:))*S2;
			
		%Remove first row/column of S1 and first row of S2(equivalent to setting E_1 to zero for now)
		S1(1,:) = [];
		S1(:,1) = [];
		S2(1,:) = [];
		
		%Combine the S1 and S2 matrices
		Sinv = S1\S2;
		
		%Save the Sinv matrix and grid parameters to the hashnamed data
		%file
		grid_params = input_params; %#ok<NASGU>
		save([datadir '/' input_hash '.mat'],'Sinv','grid_params','nodes','cubes')
		
	%If not recalculating, load Ainv
	else
		
		load([datadir '/' input_hash '.mat'],'Sinv','nodes','cubes')
		
	end
	
	
end

function [E] = processor(V,Sinv)
% Get the value of the potential function at the X,Y points

	% Turn the vector input into columns
	Vcolumns = grid_to_columns(V);
	E = [0; Sinv*Vcolumns(:)];

	E = reshape(E,size(V{1}));

end

function [gradE] = postprocessor(grid,E)
%Generate the correct outputs

	%get the dimensionality
	n_dim = numel(grid);

	% Generate the vectors for mapping ndgrid functions to gradient's
	% inputs, which expect meshgrid
	inputorder = [2 1 3:n_dim];

	
	% Prime the output
	gradE = cell(size(grid));
	gradbasis = cellfun(@(x) unique(x),grid,'UniformOutput',false);
	%Get the gradient of the scalar potential
	[gradE{inputorder}] = gradient(E,gradbasis{inputorder});

end



