function output = sysf_twin_rotor(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Twin Rotor 1-D Test'; % Display name

		case 'dependency'

			output.dependency = {};
            
		case 'initialize'
            
            %%%
            %Local Connection
			s.A_num = @(a1, a2) [zeros(size(a1)) zeros(size(a1));
                                 zeros(size(a1)) zeros(size(a1));
                                 ones(size(a1)) -ones(size(a1))];
            %%%
            % Configuration space
            s.conf_space = LieGroups.SO3;

			%%%
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1]*2.5;

			%densities for various operations
			s.density.vector = [10 10]; %density to display vector field
			s.density.scalar = [51 51]; %density to display scalar functions
			s.density.eval = [31 31];   %density for function evaluations
            s.density.finite_element=31;

			%%%
			%Display parameters

			%shape space tic locations
			s.tic_locs.x = [-2 0 2 4];
			s.tic_locs.y = [-2 0 2 4];

			%Don't optimize the reference point (turn this off for
			%non-carlike systems)
			s.xy_no_opt = 1;
            
            % Set system type variable for gait optimization
            s.system_type = 'inertia';
            
			%Output the system properties
			output = s;

	end

end