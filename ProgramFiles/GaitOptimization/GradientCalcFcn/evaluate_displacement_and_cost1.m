function [net_disp_orig, net_disp_opt, cost] = evaluate_displacement_and_cost1(s,p,tspan,ConnectionEval,IntegrationMethod,resolution)
% Evaluate the displacement and cost for the gait specified in the
% structure GAIT when executed by the system described in the structure
% S.
%
% S should be a sysplotter 's' structure loaded from a file
% sysf_FILENAME_calc.mat (note the _calc suffix)
%
% P should be a structure with fields "phi_def" and "dphi_def", returning a
% vector of shapes and shape velocities respectively. If it is not
% convenient to analytically describe the shape velocity function,
% gait.dphi should be defined as 
%
% p.dphi =  @(t) jacobianest(@(T) p.phi (T),t)
%
% as is done automatically by sysplotter, but note that this will be slower
% than specifying a function directly
%
% ConnectionEval can specify whether the local connection should be generated from
% its original function definiton, or by interpolation into the evaluated
% matrix, but is optional. Valid options are 'functional' or
% 'interpolated'. Defaults to "interpolated", which significantly faster
% when calculating the local connection or metric from scratch takes
% appreciable computational time
%
% IntegrationMethod can specify whether ODE45 or a fixed point
% (euler-exponential) integration method should be employed. Defaults to
% fixed point, to reduce interpolation overhead computational times.
%
% RESOLUTION specifies the number of points for fixed-point resolution
% evaluation. A future option may support autoconvergence, but ODE
% performance with interpolated evaluation appears to be fast enough that
% refinement of fixed-point performance is on hold.
	

	% if no ConnectionEval method is specified, default to interpolated
	if ~exist('ConnectionEval','var')
		ConnectionEval = 'interpolated';
	end
    
    % if no IntegrationMethod is specified, default to ODE
	if ~exist('IntegrationMethod','var')
		IntegrationMethod = 'fixed_step';
	end

    % if no resolution is specified, default to 100 (this only affects
    % fixed_step integration)
	if ~exist('resolution','var')
		resolution = 100;
	end

    % if p is the Fourier gait parameter, replace p to y and generate new p.
    if (~isstruct(p))&&(ismatrix(p))
        y = p;
        p = makeGait(y);
        tspan = [0 2*pi/y(end,1)];
    end

    
    
	switch IntegrationMethod
		
		case 'fixed_step'
			
			[net_disp_orig, cost] = fixed_step_integrator(s,p,tspan,ConnectionEval,resolution);
        
        case 'ODE'

            % Calculate the system motion over the gait
            sol = ode45(@(t,y) helper_function(t,y,s,p,ConnectionEval),tspan,[0 0 0 0]');

            % Extract the final motion
            disp_and_cost = deval(sol,tspan(end));

            net_disp_orig = disp_and_cost(1:3);
            cost = disp_and_cost(end);
            
        otherwise
			error('Unknown method for integrating motion');
	end

	
	% Convert the final motion into its representation in optimal
	% coordinates
    startshape = zeros(size(s.grid.eval));
	startshape_def = readGait(p.phi_def,0);
        
    actual_size = min(numel(startshape),numel(startshape_def));
    startshape(1:actual_size) = startshape_def(1:actual_size);
    
    % Interpolate beta to the start point at shape space
    startshapelist = num2cell(startshape);
    beta = cellfun(@(C) interpn(s.grid.eval{:},C,...
        startshapelist{:},'spline'),s.B_optimized.eval.Beta,...
        'UniformOutput',false);
    beta = cell2mat(beta);

    % Matrix representation of beta and the net displacement in original
    % coordinate.
    Brot = vec_to_mat_SE2(beta);
    g = vec_to_mat_SE2(net_disp_orig);

    % Adjoint action of beta to the net displacement in original coordinate
    % to acheive the net displacement in optimized coordinate.
    new_g = Brot^-1*g*Brot;
    net_disp_opt = mat_to_vec_SE2(new_g).';

	
end