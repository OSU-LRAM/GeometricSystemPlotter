function d = dhypercube(p,negcorner,poscorner)
% distance function for a hypercube, to be used by the distmesh functions

	% Get the number of dimensions
	n_dim = numel(poscorner);
	
	% Behavior based on dimensionality
	switch n_dim
		
		case 2
			
			d = drectangle0(p,negcorner(1),poscorner(1),negcorner(2),poscorner(2));
			
		case 3
			
			d = dblock0(p,negcorner(1),poscorner(1),negcorner(2),poscorner(2));
			
		otherwise
			
			error('dhypercube not yet implemented for more than 3 dimensions')
			
	end
	


	
end

