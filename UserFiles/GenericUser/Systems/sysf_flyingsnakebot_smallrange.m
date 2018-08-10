function output = sysf_flyingsnakebot_smallrange(input_mode,pathnames)

	% Default arguments
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Flying SnakeBot: 4.7 Curvature Limit'; % Display name

		case 'dependency'
%need to change these?
			output.dependency = fullfile(pathnames.sysplotterpath,...
                {'Utilities/curvature_mode_toolbox/backbone_from_curvature_bases.m',...
				'Utilities/curvature_mode_toolbox/curvatures/serpenoid_1.m',...
				'Utilities/curvature_mode_toolbox/curvatures/serpenoid_2.m',...
				'Utilities/LowRE_toolbox/LowRE_dissipation_metric_from_curvature_bases.m',...
				'Utilities/LowRE_toolbox/Inertial_local_connection_from_curvature_bases.m'});

		case 'initialize'
            
            %Functional representation of backbone shape
            s.BackboneShape = {@constant_curvature_1;@constant_curvature_2};
            
			%Functional representation of local connection
			s.A_num = @Conn_num;

            %Range of motion toggle
            rangeofmotion=4.7;
            
			%%%
			%Processing details

			%Range over which to evaluate connection
            %The 1.1 is to hide the ugly boundary vectors
			%s.grid_range = [-1,1,-1,1]*12;
            s.grid_range = [-1,1,-1,1]*(rangeofmotion+0.5);

			%densities for various operations
			s.density.vector = [11 11]; %density to display vector field
			s.density.scalar = [21 21]; %density to display scalar functions
			s.density.eval = [51 51];   %density for function evaluations   %%%%%%%%can change this for gaits 51?
			s.finite_element_density = 11;
% change this from LowRE to Inertial!!
%oh it looks like there's only "inTERNal" hmmm
%it runs while commented out still. but what are the default values?
			% power metric
			%s.metric = @(x,y) LowRE_dissipation_metric_from_curvature_bases...
			%	({@serpenoid_1;@serpenoid_2},[x;y],1,1,2);  % @(x,y) eye(2);%

			%%%
			%Display parameters

			%shape space tic locations
			%s.tic_locs.x = [-2 -1 0 1 2]*6;
			%s.tic_locs.y = [-2 -1 0 1 2]*6;
			s.tic_locs.x = [-2 -1 0 1 2]*rangeofmotion/2;
			s.tic_locs.y = [-2 -1 0 1 2]*rangeofmotion/2;


			%%%%
			%Save the system properties
			output = s;


	end

end

%This is the connection funciton evaluated at a certain point!!!
function [Ar]=Conn_num(a1,a2)

        %Original:
        % Inputs:kappa_basis_input,r,L,c,drag_ratio
		%Ar =  LowRE_local_connection_from_curvature_bases({@serpenoid_1;@serpenoid_2},[a1;a2],1,1,2);
        
        %Modified:
        % Inputs:  kappa_basis_input,r,L,density
        %OH I can only get it to run without crashing if I just don't
        %specify a density at all. I guess that works?
 		Ar =  Inertial_local_connection_from_curvature_bases({@constant_curvature_1;@constant_curvature_2},[a1;a2],1); %length of snake (should stay 1 normally)
        
end