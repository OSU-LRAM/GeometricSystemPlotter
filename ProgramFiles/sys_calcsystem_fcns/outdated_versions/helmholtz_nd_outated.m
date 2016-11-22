function [gradE, E]...
	= helmholtz_nd(grid,V,density,weight)

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
	[Sinv,nodes,tris] = preprocessor(grid,density,weight);
	
	%Multiply Ainv by the input vector field to get the scalar potential
	%function solution that has E_1 = 0 (as specified in the preprocessor)
	[E] = processor(grid,V,Sinv,nodes,tris);


	%Process E_test to get the desired outputs
	[gradE] = postprocessor(grid,E);
	
end
	

function [Sinv,nodes,tris] = preprocessor(grid,density,weight)

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
	input_params = [range(gridcolumns,1) density weight(:)'];
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
		%Construct a simplex mesh
		
		
		% Find the corner locations
		corners = find_corners(min(gridcolumns,[],1),max(gridcolumns,[],1));
		
		
		% Get the corresponding normalized mesh
		[nodes,tris] = hypercube_mesh(corners,density);

		
		% Turn the node array into cells of columns, so they can be
		% arbitrarily listed
		nodes_cell = mat2cell(nodes,size(nodes,1),ones(1,size(nodes,2)));
		
		% Interpolate the scaling
		weighti = interpn(grid{:},weight,nodes_cell{:});
		
		%Separate out useful values from grid size
		N = size(nodes,1); %number of nodes
		n_dim = size(gridcolumns,2); % number of dimensions over which fields are defined

		%Prime the S1 and S2 matrices
		S1 = zeros(N);
		S2 = zeros(N,n_dim*N);
		
		%%%%%
		%Build up the S1 and S2 matrices
		
		% Loop over each node
		for i = 1:N
			
			% Get the list of triangles containing this node
			tri_list = tris(any(tris==i,2),:);
						
			% Create a vector to hold the basis gradients for the triangles
			% as a row
			gradPhi = zeros(1,n_dim*size(tri_list,1));
						
			% Create matrices to hold the dependences of the potential
			% gradient and averaged input vector field on the values of the
			% potential function and input vector field at the nodes
			P = zeros(size(gradPhi,2),N);
			R = zeros(size(gradPhi,2),n_dim*N);
			
			% Create a matrix to hold the area of each triangle
			A = zeros(size(gradPhi,2));
			
			% Loop over each triangle
			for j = 1:size(tri_list,1)
				
				% Canonicallize the triangle by placing the ith node in the
				% first element
				tri_canon = circshift(tri_list(j,:),[0 -(find(tri_list(j,:)==i)-1)]);
							
				% Create the Jacobian from the first-quadrant, right-angle unit
				% prototype finite element to the actual triangle
				J = nodes(tri_canon(2:end),:)-repmat(nodes(tri_canon(1),:),n_dim,1);				
				
				% Fill in the basis gradient
				gradPhi_local = -ones(n_dim,1); % Sloping away from the current node, which is at (0,0)
				gradPhi((1:n_dim)+n_dim*(j-1)) = (J \ gradPhi_local)'; % Transformed from local coordinates
								
				% Fill in the dependency of the potential-function gradient
				% on the potential function at the nodes
				P_local = zeros(n_dim,N);
				for k = 1:n_dim
					P_local(k,[tri_canon(1) tri_canon(k+1)]) = [-1 1];
				end
				P((1:n_dim)+n_dim*(j-1),:) = J \ P_local;
								
				% Fill in the dependency of the averaged input vector field
				% on its value at the nodes -- simple average corresponds
				% to geometric average over the element
				for k = 1:n_dim
					R(k+n_dim*(j-1),N*(k-1) + tri_canon) = 1/length(tri_canon);
				end
				
				% Fill in the weighted area of each triangle
				for k = 1:n_dim
					A(k+n_dim*(j-1),k+n_dim*(j-1)) = det(J)/2 * mean(weighti(tri_canon)); 
				end

			end
			
			%Generate zero matrices of equal size to the permutation
			%matrices, so that the latter can be targeted at either half of
			%the input and output vectors
			Pz = zeros(size(P));
			Rz = zeros(size(R));
			
			%fill in the ith rows of S1 and S2
			S1(i,:) = gradPhi * A * P;
			S2(i,:) = gradPhi * A * R;

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
		save([datadir '/' input_hash '.mat'],'Sinv','grid_params','nodes','tris')
		
	%If not recalculating, load Ainv
	else
		
		load([datadir '/' input_hash '.mat'],'Sinv','nodes','tris')
		
	end
	
	
end

function [E] = processor(grid,V,Sinv,nodes,tris) %#ok<INUSD>
% Get the value of the potential function at the X,Y points

	% Turn the node array into cells of columns, so they can be
	% arbitrarily listed
	nodes_cell = mat2cell(nodes,size(nodes,1),ones(1,size(nodes,2)));
	
	% Interpolate the input vector field at the nodal locations
	Vi = cellfun(@(v) interpn(grid{:},v,nodes_cell{:},'cubic')...
		,V,'UniformOutput',false);	
	
	% Calculate the potential function at the nodes
	Vcolumns = grid_to_columns(Vi);
	Ei = [0;Sinv*Vcolumns(:)];
	
	% Reverse-interpolate
	n_dim = numel(grid);
	gridcolumns = grid_to_columns(grid);

	% Dilate the nodes slightly to make sure that the interpolation
	% works right
	 centerdiff = repmat((max(nodes,[],1)+min(nodes,[],1))/2,size(nodes,1),1);
	 nodes_dilated = ((nodes-centerdiff)*(1+sqrt(eps)))+centerdiff;

	if n_dim <= 3
		
		Ef = TriScatteredInterp(nodes_dilated,Ei(:),'natural');
		E = Ef(gridcolumns);
	else
		
		E = griddatan(nodes_dilated,Ei(:),gridcolumns,'linear');
		
	end
	
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



