function [nodes,tris] = hypercube_mesh(corners,density)
% Generate a normalized tetrahedral mesh over a hypercubic region with a
% given density

	% find the all-minimum corner
	mincorner = corners(end,:);

	% add this corner to the corners
	corners_normalized = corners - repmat(mincorner,size(corners,1),1);

	% find the all-maximum corner in the shifted values
	maxcorner = corners_normalized(1,:);

	% divide the corners by the largest coordinate on the maxcorner
	corners_normalized = corners_normalized/max(maxcorner);

	% hash the normalized maxcorner with the specified density
	input_params = [corners_normalized(1,:) density];
	mesh_hash = hash(input_params,'md5');

	% look for an existing copy of this mesh

	meshdir = 'sysplotter_data/mesh_data';

	%If the mesh precalc data directory doesn't exist, make it
	if ~exist(meshdir,'dir')

		mkdir(meshdir);

	end

	D3 = dir(meshdir); %list of precalculated mesh files
	
	% Get the date of this file
	thisdir = fileparts(which(mfilename));
	D4 = dir(thisdir);

	%get the file dates
	mesh_date = getfiledate(D3,[mesh_hash '.mat']);
	this_date = getfiledate(D4,[mfilename '.m']);

	%If the data file has been updated more recently than the
	%helmholtz function, don't need to rerun the calculation
	remesh = 1;
	if (mesh_date > this_date)

		%Also, check the actual size and spacing to protect against hash
		%collisions
		load([meshdir '/' mesh_hash '.mat'],'grid_params')
		params_test = sum(grid_params == input_params); %#ok<NODEF>
		if params_test == length(grid_params)

			remesh = 0;

		end

	end

	% If we need to remesh, recalculate the mesh in normalized
	% coordinates and cache it
	if remesh

		% define the n-dimensionalbox bounding the target region of the shape space, so
		% that points outside this box can be discarded
		fd = @(p) dhypercube0(p,corners_normalized(end,:),corners_normalized(1,:));

		% Find the target length of the simplex edges
		start_length = min(corners_normalized(1,:))/density;

		% Create the simplex mesh
		[nodes,tris]=distmeshnd_hypercube(fd,@huniform,start_length...
			,corners_normalized([end 1],:)...
			,corners_normalized...
			,[]);

		[nodes,tris] = fixmesh(nodes,tris);

		% Save out the simplex mesh
		grid_params = input_params; %#ok<NASGU>
		save([meshdir '/' mesh_hash '.mat'],'grid_params','nodes','tris')

	else

		load([meshdir '/' mesh_hash '.mat'],'nodes','tris')

	end
	
	%%%%%%%%%%%%%%%
	% Rescale the mesh to its original size
	nodes = nodes*max(maxcorner);
	
	% reposition the mesh
	nodes = nodes + repmat(mincorner,size(nodes,1),1);
	
end