function output = sysf_diffdrive(input_mode,pathnames)
%Return the system properties for a diffdrive car

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'name';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Diffdrive Car'; % Display name

		case 'dependency'

			output.dependency = {}; 
			% For system files depending on other files in the Systems
			% directory, list those filenames here for the dependency
			% checker

		case 'initialize'

			% Differential dr
            
			%%%
			% Local connection (functions at end of file)
			s.A_num = @Conn_num;


			%%%
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1]*2.5;

			%densities for various operations
			s.density.vector = [10 10]; %density to display vector field
			s.density.scalar = [51 51]; %density to display scalar functions
			s.density.eval = [21 21];   %density for function evaluations
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

			%%%%
            
			%Output the system properties
			output = s;
		
	end
    
end

function A_num = Conn_num(a1,a2)
% Numerator of the local connection (in this case, no denominator is needed). 
%
% Note that the function returns a single large array of values for the
% local connection:
% 
% if a1 and a2 are each NxM matrices of joint angles, A is a 3x2 collection
% of NxM matrices, each corresponding to the i,jth component of A

	A_num = [-ones(size(a1)) -ones(size(a1));
		zeros(size(a1)) zeros(size(a1));
		ones(size(a1)) -ones(size(a1))];
	
end
