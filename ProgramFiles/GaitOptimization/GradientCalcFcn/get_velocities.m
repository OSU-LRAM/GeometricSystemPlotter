% Evaluate the body velocity and cost velocity (according to system metric)
% at a given time
function [gcirc, dcost] = get_velocities(t,s,gait,ConnectionEval,A,metric,dM)

	% Get the shape and shape derivative at the current time
    shape = zeros(size(s.grid.eval));
    dshape = zeros(size(s.grid.eval));
    ddshape = zeros(size(s.grid.eval));
    
	shape_gait_def = readGait(gait.phi_def,t);
	dshape_gait_def = readGait(gait.dphi_def,t);
    ddshape_gait_def = readGait(gait.ddphi_def,t);
    
    actual_size = min(numel(shape),numel(shape_gait_def));
    shape(1:actual_size) = shape_gait_def(1:actual_size);
    dshape(1:actual_size) = dshape_gait_def(1:actual_size);
    ddshape(1:actual_size) = ddshape_gait_def(1:actual_size);
  
    M_a = metric;
    
	shapelist = num2cell(shape);
	
    % If doing functional eval of system (not recommended)
	% Get the local connection and metric at the current time, in the new coordinates
	if strcmpi(ConnectionEval,'functional')
			
        A = s.A_num(shapelist{:})./s.A_den(shapelist{:});

        switch s.system_type
            case 'drag'
                metric = s.metric(shapelist{:});
            case 'inertia'
                error('Functional ConnectionEval method not supported for inertia systems!')
        end

    end
	
	% Get the body velocity at the current time
	%t;
    gcirc = - A * dshape(:);

    switch s.costfunction
        case {'pathlength metric','pathlength coord'}
            dcost = sqrt(dshape(:)'*metric*dshape(:));
        case 'pathlength metric2'
            dcost = sqrt(dshape(:)'*metric*metric*dshape(:));
        case 'torque'
            dcost = torque_cost(M_a,dM,shape,dshape,ddshape,metric);
        case 'covariant acceleration'
            dcost = acceleration_cost(M_a,dM,shape,dshape,ddshape,metric);
        case 'acceleration coord'
            dcost = ddshape(:)'*metric*ddshape(:);
        case 'power quality'
            dcost = power_quality_cost(M_a,dM,shape,dshape,ddshape);
    end
	
end