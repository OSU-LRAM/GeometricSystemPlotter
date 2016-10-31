function [gradE, E]...
	= helmholtz_nd(grid,V,weight)

% Discrete Hodge-Helmholtz-esqe decomposition on a rectangular region, for optimal x,y point, with
% optimization to store the transform matrix. 
%
% Inputs are:
%  grid: cell of ndgridoutputs
%  Vx,Vy: vector components at the meshgrid points
%  density: number of nodes to use along shortest dimension
%  weight: relative importance of each X,Y point (optional)
%
% Outputs are: 
%  gradE: Curl-free component of the vector field
%  W: Remainder of the vector field (approximately divergence free)
%  E: Scalar potential function whose gradient is the curl-free component,
%  and which has a value of zero at the top left corner
%

	%Get the matrix mapping from velocities to projected scalar potential
	if ~exist('weight','var')
		weight = ones(size(grid{1}));
	end
	[Sinv,nodes,tris] = preprocessor(grid,weight);
	
	%Multiply Ainv by the input vector field to get the scalar potential
	%function solution that has E_1 = 0 (as specified in the preprocessor)
	[E] = processor(grid,V,Sinv,nodes,tris);


	%Process E_test to get the desired outputs
	[gradE] = postprocessor(grid,E);
	
end
	

function [Sinv,nodes,cubes] = preprocessor(grid,weight)

	% get the location of this file
	[pathstr] = fileparts(which([mfilename '.m']));
	addpath([pathstr, '/distmesh'])
	
	%%%%%%%%%
	%First, check if we've already calculated this matrix

	recalc = 1; %default value of flag for recalculation
		
	%Name the directory in which we store precalculated data
	datadir = ['sysplotter_data/' mfilename '_data'];
	
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
		
		% Loop over each node
		for i = 1:N
			
			% Get the list of cubes containing this node
			cube_list = cubes(any(cubes==i,2),:);
						
			% Create a vector to hold the basis gradients for the cubes
			% as a row
			gradPhi = zeros(1,n_dim*size(cube_list,1));
						
			% Create matrices to hold the dependences of the potential
			% gradient and averaged input vector field on the values of the
			% potential function and input vector field at the nodes
			P = zeros(size(gradPhi,2),N);
			R = zeros(size(gradPhi,2),n_dim*N);
			
			% Create a matrix to hold the area of each cube
			A = zeros(size(gradPhi,2));
			
			
			% Loop over each cube attached to the node
			for j = 1:size(cube_list,1)
				
				% Get which corner (and thus basis function) of the cube is
				% the current node
				cornerI = find(cube_list(j,:)==i);
					
				%%%%%
				% Fill in the basis gradient
				
				% basis function derivative in local coordinates for the
				% basis function that is non-zero at this node, evaluated
				% at the quadrature points
				gradPhi_local_quad = [dshape_at_quad{cornerI}{:}]';
				
				% Evaluate the jacobian for this element at the integration
				% points
				J_element = arrayfun(@(varargin) J(varargin,j),quad_points{:},'UniformOutput',false);
				
				% transform derivative into real coordinates
				gradPhi_transformed_quad = cell2mat(...
					cellfun(@(x,y) x\y, J_element',...      % transpose is of the cell array, not the contents
					mat2cell(gradPhi_local_quad,size(gradPhi_local_quad,1)...
						,ones(1,size(gradPhi_local_quad,2)) ),'UniformOutput',false) );
				
				% dot-product the quadrature weights with the gradient at
				% the quadrature points to get the integrated basis gradient
				% value for the element
				gradPhi_transformed = gradPhi_transformed_quad*quad_weights;
				
				% Place into the overall gradPhi vector for this node
				gradPhi((1:n_dim)+n_dim*(j-1)) = gradPhi_transformed;
				
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
						cellfun(@(x,y) x\y, J_element',... % transpose is of the cell array, not the contents
						mat2cell(P_local_quad{m},size(P_local_quad{m},1)...
							,ones(1,size(P_local_quad{m},2)) ),'UniformOutput',false) )';

				end

				% dot-product the quadrature weights with the gradient at
				% the quadrature points to get the integrated potential gradient
				% value for the element. 
				P_transformed = cellfun(@(p) quad_weights'*p,P_transformed_quad,'UniformOutput',false);

				% Slot the gradient values into a vector that pairs them
				% with the nodal values they correspond to
				for k = 1:length(cube_list(j,:))
					
					P((1:n_dim)+n_dim*(j-1),cube_list(j,k)) = P_transformed{k}';

				end
								
				% Fill in the dependency of the averaged input vector field
				% on its value at the nodes 
				
				% In local coordinates, this dependency is the value of the
				% shape functions at the quadrature points
				R_local_quad = cellfun(@(x)x',shape_at_quad,'UniformOutput',false);
				
				% In real coordinates, this value is scaled by the area of
				% the element, encoded as the jacobian's determinant
				R_transformed_quad = cell(size(R_local_quad));
				for m = 1:numel(R_transformed_quad)

					R_transformed_quad{m} = cell2mat(...
						cellfun(@(x,y) det(x)*y, J_element',... % transpose is of the cell array, not the contents
						mat2cell(R_local_quad{m},size(R_local_quad{m},1)...
							,ones(1,size(R_local_quad{m},2)) ),'UniformOutput',false) )';

				end
				
				% dot-product the quadrature weights with the vector dependency at
				% the quadrature points to get the integrated vector value
				% value for the element. 
				R_transformed = cellfun(@(r) quad_weights'*r,R_transformed_quad,'UniformOutput',false);
				
				
				% Slot the gradient values into a vector that pairs them
				% with the nodal values they correspond to (this is n_dim
				% slots, one for each vector component
				for k = 1:length(cube_list(j,:))
					
					for m = 1:n_dim
						R(m+n_dim*(j-1),N*(m-1)+cube_list(j,k)) = R_transformed{k};
					end

				end
				

			end
						
			%fill in the ith rows of S1 and S2
			S1(i,:) = gradPhi * P;
			S2(i,:) = gradPhi * R;

		end
			
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

function [E] = processor(grid,V,Sinv,nodes,tris) %#ok<INUSD>
% Get the value of the potential function at the X,Y points

	% Turn the vector input into columns
	Vcolumns = grid_to_columns(V);
	E = [0; Sinv*Vcolumns(:)];

	E = reshape(E,size(V{1}));
% 
% 	% Turn the node array into cells of columns, so they can be
% 	% arbitrarily listed
% 	nodes_cell = mat2cell(nodes,size(nodes,1),ones(1,size(nodes,2)));
% 	
% 	% Interpolate the input vector field at the nodal locations
% 	Vi = cellfun(@(v) interpn(grid{:},v,nodes_cell{:},'cubic')...
% 		,V,'UniformOutput',false);	
% 	
% 	% Calculate the potential function at the nodes
% 	Vcolumns = grid_to_columns(Vi);
% 	Ei = [0;Sinv*Vcolumns(:)];
% 	
% 	% Reverse-interpolate
% 	n_dim = numel(grid);
% 	gridcolumns = grid_to_columns(grid);
% 
% 	% Dilate the nodes slightly to make sure that the interpolation
% 	% works right
% 	 centerdiff = repmat((max(nodes,[],1)+min(nodes,[],1))/2,size(nodes,1),1);
% 	 nodes_dilated = ((nodes-centerdiff)*(1+sqrt(eps)))+centerdiff;
% 
% 	if n_dim <= 3
% 		
% 		Ef = TriScatteredInterp(nodes_dilated,Ei(:),'natural');
% 		E = Ef(gridcolumns);
% 	else
% 		
% 		E = griddatan(nodes_dilated,Ei(:),gridcolumns,'linear');
% 		
% 	end
% 	
% 	E = reshape(E,size(V{1}));

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



