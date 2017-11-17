function output = sysf_serpenoid_extendable_0(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
    end
		
    % curvature function things. constant, function definitions and name
    stretch_const = 0;
    stretchdef = [num2str(stretch_const),'*(a1*cos(2*pi*s)+a2*sin(2*pi*s))^2+1'];
    curvdef = 'a1*cos(2*pi*s)+a2*sin(2*pi*s)';
    % creates file name with b constant 1000 times the actual, so there's
    % no decimal issues
    curvdef_name = ['serpenoid_stretch_b_',num2str(stretch_const*1000)]; 
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = ['Serpenoid Stretchable Snake, b = ',num2str(stretch_const)]; % Display name

		case 'dependency'

			output.dependency = [fullfile(pathnames.sysplotterpath,...
                {'Utilities/curvature_mode_toolbox/backbone_from_stretchable_curvature.m',...
				'Utilities/curvature_mode_toolbox/curvatures/serpenoid_1.m',...
                'Utilities/curvature_mode_toolbox/curvatures/serpenoid_2.m',...
				'Utilities/LowRE_toolbox/LowRE_dissipation_metric_from_stretchable_curvature.m',...
				'Utilities/LowRE_toolbox/LowRE_local_connection_from_stretchable_curvature.m'}),...
                fullfile(pathnames.syspath,curvdef_name,'.mat')];

		case 'initialize'

            % generate curvdef function for stretchable serpenoid motion
            % given the b constant defined above
            
            make_curvdef_stretch({curvdef,stretchdef},{'a1','a2'},...
                ['serpenoid_stretch_b_',num2str(stretch_const*1000)],[1 1 1 1 1 1],...
                pathnames.sysplotterpath,pathnames.syspath)

            curvdef_fun = str2func(['curv_',curvdef_name]);
            
			% Functional representation of local connection
			s.A_num = @(a1,a2) Conn_num(a1,a2,curvdef_fun);


			%%%
			% Processing details

			% Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1]*20; %*2.5;
            % s.grid_range = [-1,1,-1,1]*ceil((3/0.1)*stretch_const + 15);

			% densities for various operations
			s.density.vector = [11 11]; % density to display vector field
			s.density.scalar = [21 21]; % density to display scalar functions
			s.density.eval = [21 21];   % density for function evaluations
			s.finite_element_density = 11;
			% power metric
			s.metric = @(x,y) M_helper(x,y,curvdef_fun);  % @(x,y) eye(2);%

			%%%
			%Display parameters

			% shape space tic locations
			s.tic_locs.x = [-2 -1 0 1 2]*6; %*1;
			s.tic_locs.y = [-2 -1 0 1 2]*6; %*1;


			%%%%
			% Save the system properties
			output = s;


	end

end

function [Ar]=Conn_num(a1,a2,curvdef_fun)

    Ar = LowRE_local_connection_from_stretchable_curvature(curvdef_fun,[a1;a2],1,1,2);

end

function M = M_helper(x,y,curvdef_fun)

    M = LowRE_dissipation_metric_from_stretchable_curvature(curvdef_fun,[x;y],1,1,2);

end
