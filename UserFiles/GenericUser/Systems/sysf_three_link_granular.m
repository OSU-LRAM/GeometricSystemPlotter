function output = sysf_three_link_granular(input_mode,pathnames)
% This file is an example of how to load the local connection as a data
% file, instead of generating it programmatically.


	% Default arguments
    if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
    end
	

    % File name where local connection is stored. This should be a file
    % with:
    %
    % alpha1 and alpha2: ndgrid matrices (*not* meshgrid matrices)
    % Ax1 through Atheta2: each component of the local connection,
    %       evaluated at all grid points
    Local_connection = 'Granular_Local_Connection_Matrix.mat';

    % File name where power-usage data is stored
    %Metric = 'Granular_Metric_Tensor.mat'; %Uncomment this line if loading a metric tensor

	
	%%%%%%%
	
	switch input_mode

		case 'name'

			output = 'Granular Swimmer'; % Display name

		case 'dependency'

			output.dependency = {fullfile(pathnames.syspath,Local_connection)};
            
		case 'initialize'

			%Initialize a kinematic snake with unit values for L and R

			%%%%%
			% Filename to save to
			output = mfilename;

            
            % Load the connection and metric data
            load(Local_connection)
            %load(Metric)

			%Functional representation of local connection
			s.A_num = @(a1,a2) Conn_num(a1,a2,alpha1,alpha2,Ax1,Ax2,Ay1,Ay2,Atheta1,Atheta2);

			%%%
            
			%Processing details

			%Mark that system has a singularity that needs to be accounted for
			s.singularity = 0;

			%Range over which to evaluate connection
			s.grid_range = [-1,1,-1,1]*2.5;

			%densities for various operations
			s.density.vector = [11 11]; %density to display vector field
			s.density.scalar = [21 21]; %density to display scalar functions
			s.density.eval = [11 11];   %density for function evaluations
            s.density.metric_eval = [11 11]; %density for metric evaluation
            s.density.finite_element=31;

            %%%%%%
            % Define system geometry
            s.geometry.type = 'n-link chain';
            s.geometry.linklengths = [1 1 1];
            s.geometry.baseframe = 'center';
            s.geometry.length = 1;
            
            % Define properties for visualizing the system
            
            % Make a grid of values at which to visualize the system in
            % illustrate_shapespace. (Use a cell of gridpoints along each
            % axis to use different spacings for different axes)
            s.visual.grid_spacing = [-1  0  1];


            % power metric
			%s.metric = eye(2); %@(x,y) Granular_metric_calc(x,y,Metric_Tensor_raw,alpha1,alpha2);

			%%%
            
			%Display parameters

			%shape space tic locations
			s.tic_locs.x = [-1 0 1]*2;
			s.tic_locs.y = [-1 0 1]*2;

            % Set system type variable for gait optimization
            s.system_type = 'drag';
            
			%%%%
            
			%Save the system properties
			output = s;


	end

end

function [A]=Conn_num(a1,a2,alpha1,alpha2,Ax1,Ax2,Ay1,Ay2,Atheta1,Atheta2)
% This function takes the components of the connection and arranges them
% into the internal structure used by the sysplotter back-end code.
% 
% To improve the quality of the results, natural symmetries in the geometry
% of the system can be exploited to average the connection terms that
% should be the same as each other.
%
% There are three options for doing this, depending on the kind of symmetry
% expected in the data:
%
% 3link: If the system has three links with equal first and last links
% EvenOdd: If the system shape variables are for even and odd modes
% None: Don't average the connection across any points, and display exactly
% the data in the .mat file

% Select averaging mode
 averaging_mode = '3link';
% averaging_mode = 'EvenOdd';
% averaging_mode = 'None';


A = four_way_symmetry(alpha1, alpha2, Ax1, Ax2, Ay1, Ay2, Atheta1,Atheta2, a1, a2, 1/3, averaging_mode);
        

end