function output = sysf_honey_swimmer_serpenoid(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
	
	if ~exist('datapath','var')
		
		datapath = '';
	end
	
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Low Re Serpenoid first modes'; % Display name

		case 'dependency'

			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Utilities/curvature_mode_toolbox/backbone_from_curvature_bases.m',...
				'Utilities/curvature_mode_toolbox/curvatures/serpenoid_1.m',...
				'Utilities/curvature_mode_toolbox/curvatures/serpenoid_2.m',...
				'Utilities/LowRE_toolbox/LowRE_dissipation_metric_from_curvature_bases.m',...
				'Utilities/LowRE_toolbox/LowRE_local_connection_from_curvature_bases.m'});

		case 'initialize'


			%Functional representation of local connection
			s.A_num = @Conn_num;


			%%%
			%Processing details

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1]*12;

			%densities for various operations
			s.density.vector = [11 11]; %density to display vector field
			s.density.scalar = [21 21]; %density to display scalar functions
			s.density.eval = [21 21];   %density for function evaluations
			s.finite_element_density = 11;
			% power metric
			s.metric = @(x,y) LowRE_dissipation_metric_from_curvature_bases...
				({@serpenoid_1;@serpenoid_2},[x;y],1,1,2);  % @(x,y) eye(2);%

			%%%
			%Display parameters

			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*6;
			s.tic_locs.y = [-1 0 1]*6;


			%%%%
			%Save the system properties
			output = s;


	end

end

function [Ar]=Conn_num(a1,a2)


		Ar =  LowRE_local_connection_from_curvature_bases({@serpenoid_1;@serpenoid_2},[a1;a2],1,1,2);

end