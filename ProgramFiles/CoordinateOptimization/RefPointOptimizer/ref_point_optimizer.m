function [gradE_x, gradE_y, Ex, Ey]...
	= ref_point_optimizer(grid,Vx,Vy,Vt,weight,data_storage_location)

% Discrete Reference Point Optimization for an SE(2) system. 
%
% Inputs are:
%  grid: cell of ndgridoutputs
%  V: vector field structure -- cell array of arrays, nth cell
%    corresponding to nth component
%  weight: relative importance of each point (optional, defaults to even weighting)
%  data_storage_location: directory in which to cache results of
%    calculation (optional, defaults to the current directory)
%
% Outputs are: 
%  gradE_x, GradE_y: Curl-free component of the vector fields
%  E_x, E_y: Scalar potential functions whose gradients are the curl-free component,
%  and which has a value of zero at the top left corner
%

	% Set default arguments
	if ~exist('weight','var')
		weight = ones(size(grid{1}));
	end
	
	if ~exist('data_storage_location','var')
		data_storage_location = '.';
	end
	
	%Get the matrix mapping from velocities to projected scalar potential	
	[Sinv,nodes,cubes] = preprocessor(grid,Vt,weight,data_storage_location); %#ok<NASGU,ASGLU>
	
	%Multiply Ainv by the input vector field to get the scalar potential
	%function solution that has E_1 = 0 (as specified in the preprocessor)
	[Ex, Ey] = processor(Vx,Vy,Sinv);


	%Process E_test to get the desired outputs
	[gradE_x] = postprocessor(grid,Ex);
	[gradE_y] = postprocessor(grid,Ey);
	
end
	

function [Sinv,nodes,cubes] = preprocessor(grid,Vt,weight,data_storage_location)

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
		%Construct a cuboid mesh -- assume regular underlying grid
		
		% Stages numbered according to Becker Carey Oden page 5.3.2
		
		
		
		%%%%%%%%%%%%%
		% Stage 1: Get the node list for the mesh and for each hypercube
		[nodes,cubes] = hypercube_mesh(grid);

		
		%Separate out useful values from grid size
		N = size(gridcolumns,1); %number of nodes
		n_dim = size(gridcolumns,2); % number of dimensions over which fields are defined

		%Prime the S1 and S2 matrices
		S1 = zeros(2*N);
		S2 = zeros(2*N,2*n_dim*N);
		
		% Define the (n-)linear shape functions for the finite elements
		[shape_functions, dshape_functions] = hypercube_element_shape_functions(n_dim);

		% Extract the number of shape functions
		n_shape = numel(shape_functions);
		
		%%%%%%%%%%%%%%
		% Stage 2: Get the quadrature points and weights for the elements
		[quad_points, quad_weights] = hypercube_quadrature_init(n_dim);
		
		
		%%%%%%%%%%%%
		% Stage 3: Calculate the values of the shape functions and their
		% derivatives at the quadrature points
		
		shape_at_quad = cellfun(@(f) f(quad_points),shape_functions,'UniformOutput',false);
		
		dshape_at_quad = cellfun(@(f) f(quad_points),dshape_functions,'UniformOutput',false);
		
		% Construct interpolating matrices
		shape_interpolator = cat(2,shape_at_quad{:});
		
		dshape_interpolator = cell(1,size(dshape_at_quad,2));
		for i = 1:size(dshape_interpolator,2)
			dshape_interpolator{i} = cat(2,dshape_at_quad{:,i});
		end
		
		
		%%%%%%%%%%%
		% Stage 4: Calculate the values of the problem coordinates at the
		% quadrature points. We're on a regular grid, so can just use the
		% first element
		
		quad_points_real = cell(n_dim,1);      % x(eta), etc
		quad_deriv_real = cell(n_dim,n_dim);   % dx/d eta , etc
		for i = 1:size(quad_points_real,1) % Rows are each one dimension
			
			quad_points_real{i} = shape_interpolator*nodes(cubes(1,:),i);
			
			for j = 1:size(quad_deriv_real,2) % Columns are derivatives with respect to the jth dimension
				
				quad_deriv_real{i,j} = dshape_interpolator{j}*nodes(cubes(1,:),i);
				
			end
			
		end
		
		
		%%%%%%%%%%%%%
		% Stage 5: Construct the Jacobian
		
		% Build the Jacobian at all nodes
		
		% Reshape and merge the quadrature derivatives
		quad_deriv_reshape = cell2mat(cellfun(@(x) permute(x,[3 2 1]),quad_deriv_real,'UniformOutput',false));
		
		% Take the first slice as the Jacobian, since we're on a regular
		% grid and the jacobian is thus constant.
		J = quad_deriv_reshape(:,:,1); %dx/d eta
				
		
		%%%%%%%%%%%%
		% Stage 6: Get shape function derivatives with respect to problem
		% coordinates
		
		% First, merge dshape_at_quad into a single array
		dshape_at_quad_merged = cell2mat(dshape_at_quad);
		
		% Multiply through by the jacobian
		dshape_at_quad_merged_real = (J\dshape_at_quad_merged')';
		
		% reseparate
		dshape_at_quad_real = mat2cell(dshape_at_quad_merged_real,length(quad_points{1})*ones(n_shape,1),ones(1,n_dim));
		
		% Reconstruct interpolating matrix for the derivatives
		dshape_interpolator = cell(1,size(dshape_at_quad,2));
		for i = 1:size(dshape_interpolator,2)
			dshape_interpolator{i} = cat(2,dshape_at_quad_real{:,i});
		end		

		%%%%%%%%%%%%%%%%%%%%%%%%%%%

		% Construct the elements of the solution matrices (_q is for
		% versions with quadrature components broken out

		% Quadrature integration vector
		q_dim = repmat(quad_weights',1,n_dim);

		% Del operator
		Del_q = cat(1,dshape_interpolator{:});
		
		% gradPhi operators (same as columns of Del, but splitting the name
		% for clarity below
		gradPhi_q = Del_q;
		
		% gradPhi_dot_Del operators
		gradPhi_dot_Del_q = cell(1,size(Del_q,2));
		for i = 1:size(gradPhi_dot_Del_q,2)
			gradPhi_dot_Del_q{i} = repmat(gradPhi_q(:,i),1,size(Del_q,2)).*Del_q;
		end
		% Apply quadrature integration
		gradPhi_dot_Del = cellfun(@(x) q_dim*x,gradPhi_dot_Del_q,'UniformOutput',false);
		
		% "Selection" operator -- get values of function at quad points
		Sel_q = shape_interpolator;
		Sel_q_dim = repmat(Sel_q,n_dim,1);
		
		% Phi operator (same as selection operator, but splitting for
		% clarity
		Phi_q = Sel_q;
		Phi_q_dim = Sel_q_dim;
		
		% gradPhi times a selection matrix
		gradPhi_dot_R = cell(n_dim,size(Del_q,2)); % one row per dimension, one column per shape function
		for i = 1:size(gradPhi_dot_R,1)
			for j = 1:size(gradPhi_dot_R,2)
				gradPhi_dot_R{i,j} = quad_weights(:)' *...
					(repmat(dshape_interpolator{i}(:,j),1,n_shape).* Sel_q);
			end
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
			
				% Build matrix contributions that depend on the cell we're
				% in
				
				% Vector values at nodal points
				W_selected = vecentries(cube_list(j,:),:);% cell2mat(cellfun(@(x) x(cube_list(j,:)),Vt(:),'UniformOutput',false));

				% Get the vector-field values at the quadrature points
				W_sel_q = Sel_q * W_selected;
				
				% Get the dot-product of W with itself
				W_dot_W_q = W_sel_q(:).^2;
				
				% Get Phi_W_dot_W
				phiW_dot_W_q = Phi_q_dim(:,cornerI).*W_dot_W_q;
				
				% Multiply this by the quadrature weights and another
				% selection function to generate the map from potential
				% function to matrix element
				phiW_dot_W = q_dim * (repmat(phiW_dot_W_q,1,n_shape) ...
					.* Sel_q_dim);
				
				
				% gradPhi dot W
				gradPhi_dot_W_q = gradPhi_q(:,cornerI) .* W_sel_q(:);
				gradPhi_dot_W = q_dim * (repmat(gradPhi_dot_W_q,1,n_shape) ...
					.* Sel_q_dim);
				
				
				
				% phiW_dot_del
				phiW_dot_Del_q = repmat(Phi_q_dim(:,cornerI),1,n_shape).*repmat(W_sel_q(:),1,n_shape).*Del_q;
				phiW_dot_Del = q_dim * phiW_dot_Del_q;
				
				
				% phiW_dot_R
				phiW_dot_R = cell(n_dim,1);
				for k = 1:size(phiW_dot_R,1)
					phiW_dot_R{k} = quad_weights(:)' * (repmat(Phi_q(:,cornerI),1,n_shape).*(repmat(W_sel_q(:,k),1,n_shape).* Sel_q));
				end

				
				

				% Slot the dotproduct values into a vector that pairs them
				% with the nodal values they correspond to

				% pull out current nodes for easier readability
				cn = cube_list(j,:);


				% Left hand matrix

				S1(i,	cn)		= S1(i,		cn)		+ ( gradPhi_dot_Del{cornerI}	+ phiW_dot_W); % Top, left block:
				S1(i,	cn+N)	= S1(i,		cn+N)	+ (-gradPhi_dot_W				+ phiW_dot_Del); % Top, right block

				S1(i+N,	cn)		= S1(i+N,	cn)		+ ( gradPhi_dot_W				- phiW_dot_Del); % Bottom, left block
				S1(i+N,	cn+N)	= S1(i+N,	cn+N)	+ ( gradPhi_dot_Del{cornerI}	+ phiW_dot_W); % Bottom, right block




				% Offset each component's-worth of R values by the number
				% of nodes
				for m = 1:n_dim
					S2(i,	cn+(m-1)*N)				= S2(i,		cn+(m-1)*N)				- gradPhi_dot_R{m,cornerI}; % Top left block
					S2(i,	cn+(m-1)*N+(N*n_dim))	= S2(i,		cn+(m-1)*N+(N*n_dim))	- phiW_dot_R{m}; % Top right block

					S2(i+N,	cn+(m-1)*N)				= S2(i+N,	cn+(m-1)*N)				+ phiW_dot_R{m}; % Bottom left block
					S2(i+N,	cn+(m-1)*N+(N*n_dim))	= S2(i+N,	cn+(m-1)*N+(N*n_dim))	- gradPhi_dot_R{m,cornerI}; % Bottom right block
				end


			end
					
				
				
		end
		
		% Multiply the matrices by the node weighting factors
		%S1 = diag([weight(:);weight(:)])*S1;
		S2 = diag([weight(:);weight(:)])*S2;		
        
        %Combine the S1 and S2 matrices
		Sinv = S1\S2;
		
		%Save the Sinv matrix and grid parameters to the hashnamed data
		%file
		grid_params = input_params; %#ok<NASGU>
		save(fullfile(datadir,[input_hash '.mat']),'Sinv','grid_params','nodes','cubes')
		
	%If not recalculating, load Ainv
    else
		
        try
            load(fullfile(datadir, [input_hash '.mat']),'Sinv','nodes','cubes')
            %Sometimes Sinv doesn't load if the file is too large.  This is to make it trigger an
            %error to recalc if it doesn't
            Sinv;
        catch


            %%%%%%%%
            %Construct a cuboid mesh -- assume regular underlying grid

            % Stages numbered according to Becker Carey Oden page 5.3.2



            %%%%%%%%%%%%%
            % Stage 1: Get the node list for the mesh and for each hypercube
            [nodes,cubes] = hypercube_mesh(grid);


            %Separate out useful values from grid size
            N = size(gridcolumns,1); %number of nodes
            n_dim = size(gridcolumns,2); % number of dimensions over which fields are defined

            %Prime the S1 and S2 matrices
            S1 = zeros(2*N);
            S2 = zeros(2*N,2*n_dim*N);

            % Define the (n-)linear shape functions for the finite elements
            [shape_functions, dshape_functions] = hypercube_element_shape_functions(n_dim);

            % Extract the number of shape functions
            n_shape = numel(shape_functions);

            %%%%%%%%%%%%%%
            % Stage 2: Get the quadrature points and weights for the elements
            [quad_points, quad_weights] = hypercube_quadrature_init(n_dim);


            %%%%%%%%%%%%
            % Stage 3: Calculate the values of the shape functions and their
            % derivatives at the quadrature points

            shape_at_quad = cellfun(@(f) f(quad_points),shape_functions,'UniformOutput',false);

            dshape_at_quad = cellfun(@(f) f(quad_points),dshape_functions,'UniformOutput',false);

            % Construct interpolating matrices
            shape_interpolator = cat(2,shape_at_quad{:});

            dshape_interpolator = cell(1,size(dshape_at_quad,2));
            for i = 1:size(dshape_interpolator,2)
                dshape_interpolator{i} = cat(2,dshape_at_quad{:,i});
            end


            %%%%%%%%%%%
            % Stage 4: Calculate the values of the problem coordinates at the
            % quadrature points. We're on a regular grid, so can just use the
            % first element

            quad_points_real = cell(n_dim,1);      % x(eta), etc
            quad_deriv_real = cell(n_dim,n_dim);   % dx/d eta , etc
            for i = 1:size(quad_points_real,1) % Rows are each one dimension

                quad_points_real{i} = shape_interpolator*nodes(cubes(1,:),i);

                for j = 1:size(quad_deriv_real,2) % Columns are derivatives with respect to the jth dimension

                    quad_deriv_real{i,j} = dshape_interpolator{j}*nodes(cubes(1,:),i);

                end

            end


            %%%%%%%%%%%%%
            % Stage 5: Construct the Jacobian

            % Build the Jacobian at all nodes

            % Reshape and merge the quadrature derivatives
            quad_deriv_reshape = cell2mat(cellfun(@(x) permute(x,[3 2 1]),quad_deriv_real,'UniformOutput',false));

            % Take the first slice as the Jacobian, since we're on a regular
            % grid and the jacobian is thus constant.
            J = quad_deriv_reshape(:,:,1); %dx/d eta


            %%%%%%%%%%%%
            % Stage 6: Get shape function derivatives with respect to problem
            % coordinates

            % First, merge dshape_at_quad into a single array
            dshape_at_quad_merged = cell2mat(dshape_at_quad);

            % Multiply through by the jacobian
            dshape_at_quad_merged_real = (J\dshape_at_quad_merged')';

            % reseparate
            dshape_at_quad_real = mat2cell(dshape_at_quad_merged_real,length(quad_points{1})*ones(n_shape,1),ones(1,n_dim));

            % Reconstruct interpolating matrix for the derivatives
            dshape_interpolator = cell(1,size(dshape_at_quad,2));
            for i = 1:size(dshape_interpolator,2)
                dshape_interpolator{i} = cat(2,dshape_at_quad_real{:,i});
            end		

            %%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Construct the elements of the solution matrices (_q is for
            % versions with quadrature components broken out

            % Quadrature integration vector
            q_dim = repmat(quad_weights',1,n_dim);

            % Del operator
            Del_q = cat(1,dshape_interpolator{:});

            % gradPhi operators (same as columns of Del, but splitting the name
            % for clarity below
            gradPhi_q = Del_q;

            % gradPhi_dot_Del operators
            gradPhi_dot_Del_q = cell(1,size(Del_q,2));
            for i = 1:size(gradPhi_dot_Del_q,2)
                gradPhi_dot_Del_q{i} = repmat(gradPhi_q(:,i),1,size(Del_q,2)).*Del_q;
            end
            % Apply quadrature integration
            gradPhi_dot_Del = cellfun(@(x) q_dim*x,gradPhi_dot_Del_q,'UniformOutput',false);

            % "Selection" operator -- get values of function at quad points
            Sel_q = shape_interpolator;
            Sel_q_dim = repmat(Sel_q,n_dim,1);

            % Phi operator (same as selection operator, but splitting for
            % clarity
            Phi_q = Sel_q;
            Phi_q_dim = Sel_q_dim;

            % gradPhi times a selection matrix
            gradPhi_dot_R = cell(n_dim,size(Del_q,2)); % one row per dimension, one column per shape function
            for i = 1:size(gradPhi_dot_R,1)
                for j = 1:size(gradPhi_dot_R,2)
                    gradPhi_dot_R{i,j} = quad_weights(:)' *...
                        (repmat(dshape_interpolator{i}(:,j),1,n_shape).* Sel_q);
                end
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

                    % Build matrix contributions that depend on the cell we're
                    % in

                    % Vector values at nodal points
                    W_selected = vecentries(cube_list(j,:),:);% cell2mat(cellfun(@(x) x(cube_list(j,:)),Vt(:),'UniformOutput',false));

                    % Get the vector-field values at the quadrature points
                    W_sel_q = Sel_q * W_selected;

                    % Get the dot-product of W with itself
                    W_dot_W_q = W_sel_q(:).^2;

                    % Get Phi_W_dot_W
                    phiW_dot_W_q = Phi_q_dim(:,cornerI).*W_dot_W_q;

                    % Multiply this by the quadrature weights and another
                    % selection function to generate the map from potential
                    % function to matrix element
                    phiW_dot_W = q_dim * (repmat(phiW_dot_W_q,1,n_shape) ...
                        .* Sel_q_dim);


                    % gradPhi dot W
                    gradPhi_dot_W_q = gradPhi_q(:,cornerI) .* W_sel_q(:);
                    gradPhi_dot_W = q_dim * (repmat(gradPhi_dot_W_q,1,n_shape) ...
                        .* Sel_q_dim);



                    % phiW_dot_del
                    phiW_dot_Del_q = repmat(Phi_q_dim(:,cornerI),1,n_shape).*repmat(W_sel_q(:),1,n_shape).*Del_q;
                    phiW_dot_Del = q_dim * phiW_dot_Del_q;


                    % phiW_dot_R
                    phiW_dot_R = cell(n_dim,1);
                    for k = 1:size(phiW_dot_R,1)
                        phiW_dot_R{k} = quad_weights(:)' * (repmat(Phi_q(:,cornerI),1,n_shape).*(repmat(W_sel_q(:,k),1,n_shape).* Sel_q));
                    end




                    % Slot the dotproduct values into a vector that pairs them
                    % with the nodal values they correspond to

                    % pull out current nodes for easier readability
                    cn = cube_list(j,:);


                    % Left hand matrix

                    S1(i,	cn)		= S1(i,		cn)		+ ( gradPhi_dot_Del{cornerI}	+ phiW_dot_W); % Top, left block:
                    S1(i,	cn+N)	= S1(i,		cn+N)	+ (-gradPhi_dot_W				+ phiW_dot_Del); % Top, right block

                    S1(i+N,	cn)		= S1(i+N,	cn)		+ ( gradPhi_dot_W				- phiW_dot_Del); % Bottom, left block
                    S1(i+N,	cn+N)	= S1(i+N,	cn+N)	+ ( gradPhi_dot_Del{cornerI}	+ phiW_dot_W); % Bottom, right block




                    % Offset each component's-worth of R values by the number
                    % of nodes
                    for m = 1:n_dim
                        S2(i,	cn+(m-1)*N)				= S2(i,		cn+(m-1)*N)				- gradPhi_dot_R{m,cornerI}; % Top left block
                        S2(i,	cn+(m-1)*N+(N*n_dim))	= S2(i,		cn+(m-1)*N+(N*n_dim))	- phiW_dot_R{m}; % Top right block

                        S2(i+N,	cn+(m-1)*N)				= S2(i+N,	cn+(m-1)*N)				+ phiW_dot_R{m}; % Bottom left block
                        S2(i+N,	cn+(m-1)*N+(N*n_dim))	= S2(i+N,	cn+(m-1)*N+(N*n_dim))	- gradPhi_dot_R{m,cornerI}; % Bottom right block
                    end


                end



            end

            % Multiply the matrices by the node weighting factors
            %S1 = diag([weight(:);weight(:)])*S1;
            S2 = diag([weight(:);weight(:)])*S2;		

            %Combine the S1 and S2 matrices
            Sinv = S1\S2;

            %Save the Sinv matrix and grid parameters to the hashnamed data
            %file
            grid_params = input_params; %#ok<NASGU>
            save(fullfile(datadir,[input_hash '.mat']),'Sinv','grid_params','nodes','cubes')
        end
		
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
