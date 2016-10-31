function [gradE_x, gradE_y, Ex, Ey]...
	= ref_point_optimizer_nd_first_pass_at_reghypercubegrid(grid,Vx,Vy,Vt,weight,matrix_algorithm)

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
		weight = ones(size(X));
	end
	if ~exist('matrix_algorithm','var')
		matrix_algorithm = 'even grid optimized';
	end
	[Sinv,nodes,cubes] = preprocessor(grid,Vt,weight,matrix_algorithm); %#ok<NASGU,ASGLU>
	
	%Multiply Ainv by the input vector field to get the scalar potential
	%function solution that has E_1 = 0 (as specified in the preprocessor)
	[Ex, Ey] = processor(Vx,Vy,Sinv);


	%Process E_test to get the desired outputs
	[gradE_x] = postprocessor(grid,Ex);
	[gradE_y] = postprocessor(grid,Ey);
	
end
	

function [Sinv,nodes,cubes] = preprocessor(grid,Vt,weight,matrix_algorithm)

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
	vecentries = grid_to_columns(Vt);
	gridcolumns = grid_to_columns(grid);
	input_params = [vecentries(:)' range(gridcolumns,1) weight(:)'];
	input_hash = num2str(hash(input_params,'MD5'));

	%Compare the date of the hashnamed data file with the date of this
	%function (getfiledate will do the right thing if the file doesn't
	%exist yet).
	D1 = dir(pathstr); %list of files in this directory
	D2 = dir(datadir); %list of precalculated files
	
	%get the file dates
	ref_point_date = getfiledate(D1,[mfilename '.m']);
	data_date = getfiledate(D2,[input_hash '.mat']);
	
	%If the data file has been updated more recently than the
	%refpoint function, don't need to rerun the calculation
	if (data_date > ref_point_date)
		
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
		S1 = zeros(2*N);
		S2 = zeros(2*N,2*n_dim*N);
		
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
		
		% Select full or optimized matrix construction
		switch matrix_algorithm
		
			% Full matrix construction makes no assumptions about grid
			% regularity
			case 'full'
			

				
			case 'even grid optimized'
				
				% Jacobian is the same for all grid elements and constant
				% over the cell
				J_element = J(num2cell(zeros(n_dim,1)),1);
				detJ_element = det(J_element);
				
				% basis function derivative in local coordinates for each
				% basis function
				gradPhi_local_quad = arrayfun(@(x) [dshape_at_quad{x}{:}]',1:size(cubes,2),'UniformOutput',false);

				% transform derivative into real coordinates
				gradPhi_quad = cell(size(gradPhi_local_quad));
				for i = 1:numel(gradPhi_local_quad)
					gradPhi_quad{i} = cell2mat(...
						cellfun(@(x) J_element\x,...      
						mat2cell(gradPhi_local_quad{i},size(gradPhi_local_quad{i},1)...
							,ones(1,size(gradPhi_local_quad{i},2)) ),'UniformOutput',false) );
				end

				% Fill in the basis function
				Phi_quad = arrayfun(@(x) shape_at_quad{x}',1:size(cubes,2),'UniformOutput',false);

				
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

				%%%%%%
				% Components of S1

				% full dot product basis gradient and potential
				% function del operator
				gradPhi_dot_P_quad = cell(size(gradPhi_local_quad));
				for i = 1:numel(gradPhi_local_quad)
					gradPhi_dot_P_quad{i} ...
						= cellfun(@(p) sum(gradPhi_transformed_quad{i}.*p,1)...
							,P_transformed_quad,'UniformOutput',false);
				end
				
				%%%%%
				% Components of S2

				% keep components separate for the input vector field
				gradPhi_dot_R_quad = cell(size(gradPhi_local_quad));
				for i = 1:numel(gradPhi_local_quad) 
					gradPhi_dot_R_quad{i} ...
						= cellfun(@(r) gradPhi_transformed_quad{i}.*repmat(r,[n_dim,1])...
							,R_quad,'UniformOutput',false);
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

						%%%%
						% Theta vector field components at the quadrature
						% points
						
						% Vector values at nodal points
						W_selected = cell2mat(cellfun(@(x) x(cube_list(j,:)),Vt(:),'UniformOutput',false));
						
						
						W_quad = cellfun(@(x,y) repmat(x,1,size(y,2)).*repmat(y,size(x,1),1)...
							,num2cell(W_selected,[1 size(W_selected,2)])...
							,R_quad,'UniformOutput',false);

						W_dot_W_quad = cellfun(@(x,y) x'*x*y...
							,num2cell(W_selected,[1 size(W_selected,2)])...
							,R_quad,'UniformOutput',false);
						
						
						% Dot product of basis gradient and
						% Theta vector field
						gradPhi_dot_W_quad ...
							= cellfun(@(w) sum(gradPhi_transformed_quad{cornerI}.*...
							(w),1)...
								,W_quad,'UniformOutput',false);

						% Dot product of basis-weighted Theta vector field
						% and potential function del operator
% 						phiW_dot_P_quad ...
% 							= cellfun(@(w,p) sum(w.*p,1).*Phi_quad{cornerI}...
% 								,W_quad,P_transformed_quad,'UniformOutput',false);
						W_dot_P_quad = cellfun(@(w,p) w'*p ...
							,num2cell(W_selected,[1 size(W_selected,2)])...
							,P_transformed_quad,'UniformOutput',false);
						phiW_dot_P_quad = ...
							cellfun(@(w) w.*Phi_quad{cornerI}...
							,W_dot_P_quad,'UniformOutput',false);
							
							
							
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 						% Dot product of Theta vector field and
% 						% basis-weighted theta vector field
% 						phiW_dot_W_quad ...
% 							= cellfun(@(w) sum((w.*repmat(Phi_quad{cornerI},[n_dim,1]))...
% 							.*w)...
% 								,W_quad,'UniformOutput',false);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
						% Dot product of Theta vector field and
						% basis-weighted theta vector field
						phiW_dot_W_quad ...
							= cellfun(@(w) w.*Phi_quad{cornerI}...
								,W_dot_W_quad,'UniformOutput',false);





						% Dot product of basis-weighted Theta vector field
						% with basis
						% function
						phiW_dot_R_quad ...
							= cellfun(@(w) w.*repmat(Phi_quad{cornerI},[n_dim,1])...
								,W_quad,'UniformOutput',false);

						%%%%%%%%%%%%%%%
						% Sum the dot products over the quadrature weights to get
						% the integral value's relationship to the node values in
						% this cell
						
						% S1 components
						gradPhi_dot_P = cellfun(@(x) detJ_element*x*quad_weights(:),gradPhi_dot_P_quad{cornerI});
						gradPhi_dot_W = cellfun(@(x) detJ_element*x*quad_weights(:),gradPhi_dot_W_quad);
						phiW_dot_P = cellfun(@(x) detJ_element*x*quad_weights(:),phiW_dot_P_quad);
						phiW_dot_W = cellfun(@(x) detJ_element*x*quad_weights(:),phiW_dot_W_quad);
						
						% S2 components
						gradPhi_dot_R = cell2mat(cellfun(@(x) detJ_element*x*quad_weights(:),gradPhi_dot_R_quad{cornerI},'UniformOutput',false));
						phiW_dot_R = cell2mat(cellfun(@(x) detJ_element*x*quad_weights(:),phiW_dot_R_quad,'UniformOutput',false));

						% Slot the dotproduct values into a vector that pairs them
						% with the nodal values they correspond to
						
						% pull out current nodes for easier readability
						cn = cube_list(j,:);
						
						
						% Left hand matrix
						
						S1(i,	cn)		= S1(i,		cn)		+ ( gradPhi_dot_P	+ phiW_dot_W); % Top, left block:
						S1(i,	cn+N)	= S1(i,		cn+N)	+ (-gradPhi_dot_W	+ phiW_dot_P); % Top, right block
						
						S1(i+N,	cn)		= S1(i+N,	cn)		+ ( gradPhi_dot_W	- phiW_dot_P); % Bottom, left block
						S1(i+N,	cn+N)	= S1(i+N,	cn+N)	+ ( gradPhi_dot_P	+ phiW_dot_W); % Bottom, right block
						
						
						

						% Offset each component's-worth of R values by the number
						% of nodes
						for m = 1:n_dim
							S2(i,	cn+(m-1)*N)				= S2(i,		cn+(m-1)*N)				- gradPhi_dot_R(m,:); % Top left block
							S2(i,	cn+(m-1)*N+(N*n_dim))	= S2(i,		cn+(m-1)*N+(N*n_dim))	- phiW_dot_R(m,:); % Top right block
							
							S2(i+N,	cn+(m-1)*N)				= S2(i+N,	cn+(m-1)*N)				+ phiW_dot_R(m,:); % Bottom left block
							S2(i+N,	cn+(m-1)*N+(N*n_dim))	= S2(i+N,	cn+(m-1)*N+(N*n_dim))	- gradPhi_dot_R(m,:); % Bottom right block
						end


					end
					
				end
				
			otherwise
				
				error('Unknown matrix construction algorithm')
				
		end
		
		% Multiply the matrices by the node weighting factors
		S1 = diag([weight(:);weight(:)])*S1;
		S2 = diag([weight(:);weight(:)])*S2;		%Combine the S1 and S2 matrices
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

function [Ex, Ey] = processor(Vx,Vy,Sinv) 
% Get the value of the potential function at the X,Y points

	% Turn vector fields into one long vector array
	Vx_cols = cellfun(@(x) x(:),Vx,'UniformOutput',false);
	Vy_cols = cellfun(@(x) x(:),Vy,'UniformOutput',false);
	
	V_cols = cat(1,Vx_cols{:},Vy_cols{:});
	

	% Multiply the input vector fields by the Sinv matrix to get the total
	% potential function
	E_tot = Sinv*V_cols;
	
	% Break potential function into two pieces
	Ex = reshape(E_tot(1:length(E_tot)/2),size(Vx{1}));
	Ey = reshape(E_tot((length(E_tot)/2)+1:end),size(Vx{1}));
	

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


% 				% Loop over each node
% 				for i = 1:N
% 
% 					% Get the list of cubes containing this node
% 					cubes_at_node = find(any(cubes==i,2));
% 					cube_list = cubes(cubes_at_node,:);
% 
% 
% 
% 					% Loop over each cube attached to the node
% 					for j = 1:size(cube_list,1)
% 
% 						% Get which corner (and thus basis function) of the cube is
% 						% the current node
% 						cornerI = cube_list(j,:)==i;
% 						
% 						
% 						%%%%%
% 						% Fill in the basis gradient
% 
% 						% basis function derivative in local coordinates for the
% 						% basis function that is non-zero at this node, evaluated
% 						% at the quadrature points
% 						gradPhi_local_quad = [dshape_at_quad{cornerI}{:}]';
% 
% 						% Evaluate the jacobian for this element at the integration
% 						% points
% 						J_element = arrayfun(@(varargin) J(varargin,cubes_at_node(j)),quad_points{:},'UniformOutput',false);
% 						detJ_element = cellfun(@(x) det(x),J_element,'UniformOutput',false);
% 
% 						% transform derivative into real coordinates
% 						gradPhi_transformed_quad = cell2mat(...
% 							cellfun(@(x,y) x\y, J_element',...      % transpose is of the cell array, not the contents
% 							mat2cell(gradPhi_local_quad,size(gradPhi_local_quad,1)...
% 								,ones(1,size(gradPhi_local_quad,2)) ),'UniformOutput',false) );
% 
% 						% Fill in the basis function
% 						Phi_quad = shape_at_quad{cornerI}';
% 												
% 						
% 						%%%%%				
% 						% Fill in the dependency of the potential-function gradient
% 						% on the potential function at the nodes
% 
% 						% Values to multiply by the potential function values
% 						% to get its gradient
% 						P_local_quad = arrayfun(@(c) [dshape_at_quad{c}{:}]',1:n_shape,'UniformOutput',false);
% 
% 
% 						% Transform the potential gradient into real
% 						% coordinates
% 						P_transformed_quad = cell(size(P_local_quad));
% 						for m = 1:numel(P_transformed_quad)
% 
% 							P_transformed_quad{m} = cell2mat(...
% 								cellfun(@(x,y) x\y, J_element',... % transpose is of the cell array, not the contents
% 								mat2cell(P_local_quad{m},size(P_local_quad{m},1)...
% 									,ones(1,size(P_local_quad{m},2)) ),'UniformOutput',false) );
% 
% 						end
% 						
% 
% 
% 						% The vector fields are independent of the spacing, so they
% 						% can be interpolated using simply the shape function
% 						% values at the quadrature points
% 						R_quad = cellfun(@(x)x',shape_at_quad','UniformOutput',false);
% 
% 						%%%%
% 						% Theta vector field components at the quadrature
% 						% points
% 						
% 						% Vector values at nodal points
% 						W_selected = cell2mat(cellfun(@(x) x(cube_list(j,:)),Vt(:),'UniformOutput',false));
% 						
% 						
% 						W_quad = cellfun(@(x,y) repmat(x,1,size(y,2)).*repmat(y,size(x,1),1)...
% 							,num2cell(W_selected,[1 size(W_selected,2)])...
% 							,R_quad,'UniformOutput',false);
% 
% 
% 						%%%%%%%%%%%%%%%%%%%%%%%
% 						% Get the dot product of the basis gradient with the
% 						% potential gradient and the input vector field
% 						
% 						%%%%%%
% 						% Components of S1
% 						
% 						% full dot product basis gradient and potential
% 						% function del operator
% 						gradPhi_dot_P_quad ...
% 							= cellfun(@(p) sum(gradPhi_transformed_quad.*p,1)...
% 								,P_transformed_quad,'UniformOutput',false);
% 
% 							
% 						% Dot product of basis gradient and basis-weighted
% 						% Theta vector field
% 						gradPhi_dot_phiW_quad ...
% 							= cellfun(@(w) sum(gradPhi_transformed_quad.*...
% 							(w.*repmat(Phi_quad,[n_dim,1])),1)...
% 								,W_quad,'UniformOutput',false);
% 							
% 						% Dot product of basis-weighted Theta vector field
% 						% and potential function del operator
% 						phiW_dot_P_quad ...
% 							= cellfun(@(w,p) sum((w.*repmat(Phi_quad,[n_dim,1]))...
% 							.*p,1)...
% 								,W_quad,P_transformed_quad,'UniformOutput',false);
% 							
% 						% Dot product of Theta vector field and
% 						% basis-weighted theta vector field
% 						phiW_dot_W_quad ...
% 							= cellfun(@(w) sum((w.*repmat(Phi_quad,[n_dim,1]))...
% 							.*w)...
% 								,W_quad,'UniformOutput',false);
% 							
% 							
% 						%%%%%%%%
% 						% Components of S2
% 						
% 						% Dot product of basis gradient with interpolating
% 						% function
% 						gradPhi_dot_R_quad ...
% 							= cellfun(@(r) gradPhi_transformed_quad.*repmat(r,[n_dim,1])...
% 								,R_quad,'UniformOutput',false);
% 
% 						% Dot product of basis-weighted Theta vector field
% 						% with basis
% 						% function
% 						phiW_dot_R_quad ...
% 							= cellfun(@(w) w.*repmat(Phi_quad,[n_dim,1])...
% 								,W_quad,'UniformOutput',false);
% 
% 						
% 						%%%%%%%%%%%%%%%
% 						% Sum the dot products over the quadrature weights to get
% 						% the integral value's relationship to the node values in
% 						% this cell
% 						
% 						% S1 components
% 						gradPhi_dot_P = cellfun(@(x,y) x*y*quad_weights(:),detJ_element,gradPhi_dot_P_quad);
% 						gradPhi_dot_phiW = cellfun(@(x,y) x*y*quad_weights(:),detJ_element,gradPhi_dot_phiW_quad);
% 						phiW_dot_P = cellfun(@(x,y) x*y*quad_weights(:),detJ_element,phiW_dot_P_quad);
% 						phiW_dot_W = cellfun(@(x,y) x*y*quad_weights(:),detJ_element,phiW_dot_W_quad);
% 						
% 						% S2 components
% 						gradPhi_dot_R = cell2mat(cellfun(@(x,y) x*y*quad_weights(:),detJ_element,gradPhi_dot_R_quad,'UniformOutput',false));
% 						phiW_dot_R = cell2mat(cellfun(@(x,y) x*y*quad_weights(:),detJ_element,phiW_dot_R_quad,'UniformOutput',false));
% 
% 						% Slot the dotproduct values into a vector that pairs them
% 						% with the nodal values they correspond to
% 						
% 						% pull out current nodes for easier readability
% 						cn = cube_list(j,:);
% 						
% 						S1(i,	cn)		= S1(i,		cn)		+ ( gradPhi_dot_P		+ phiW_dot_W); % Top, left block
% 						S1(i+N,	cn+N)	= S1(i+N,	cn+N)	+ ( gradPhi_dot_P		+ phiW_dot_W); % Bottom, right block
% 						
% 						S1(i+N,	cn)		= S1(i+N,	cn)		+ ( gradPhi_dot_phiW	- phiW_dot_P); % Bottom, left block
% 						S1(i,	cn+N)	= S1(i,		cn+N)	+ (-gradPhi_dot_phiW	+ phiW_dot_P); % Top, right block
% 
% 						% Offset each component's-worth of R values by the number
% 						% of nodes
% 						for m = 1:n_dim
% 							S2(i,	cn+(m-1)*N)				= S2(i,		cn+(m-1)*N)				- gradPhi_dot_R(m,:); % Top left block
% 							S2(i,	cn+(m-1)*N+(N*n_dim))	= S2(i,		cn+(m-1)*N+(N*n_dim))	- phiW_dot_R(m,:); % Top right block
% 							S2(i+N,	cn+(m-1)*N)				= S2(i+N,	cn+(m-1)*N)				+ phiW_dot_R(m,:); % Bottom left block
% 							S2(i+N,	cn+(m-1)*N+(N*n_dim))	= S2(i+N,	cn+(m-1)*N+(N*n_dim))	- gradPhi_dot_R(m,:); % Bottom right block
% 						end
% 
% 
% 					end
% 
% 				end
