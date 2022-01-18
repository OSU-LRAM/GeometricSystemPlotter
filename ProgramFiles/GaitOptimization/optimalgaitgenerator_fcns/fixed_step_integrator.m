% Function to integrate up system velocities using a fixed-step method
function [net_disp_orig, cost] = fixed_step_integrator(s,gait,tspan,ConnectionEval,resolution)

	% Duplicate 'resolution' to 'res' if it is a number, or place res at a
	% starting resolution if an automatic convergence method is selected
	% (automatic convergence not yet enabled)
	default_res = 100;
	if isnumeric(resolution)
		res = resolution;
	elseif ischar(resolution) && strcmp(resolution,'autoconverge')
		res = default_res;
	else
		error('Unexpected value for resolution');
	end
	
	% Generate the fixed points from the time span and resolution
	tpoints = linspace(tspan(1),tspan(2),res);
	tsteps = gradient(tpoints);
    
    %Prep interpolation inputs for velocities function
    shape = zeros(size(s.grid.eval));
    shape_gait_def_0 = readGait(gait.phi_def,0);
    actual_size = min(numel(shape),numel(shape_gait_def_0));
    
    samplePoints = {};
    for dim = 1:actual_size
        samplePoints{dim} = [];
    end
    
    for time = tpoints
        shape_gait_def = readGait(gait.phi_def,time);
        shape(1:actual_size) = shape_gait_def(1:actual_size);
        for dim = 1:numel(shape)
            samplePoints{dim}(end+1) = shape(dim); 
        end
    end
    
    indexList = 1:numel(tpoints);
    id = eye(actual_size);
    
    As = cellfun(@(C) -interpn(s.grid.eval{:},C,...
        samplePoints{:},'spline'),s.vecfield.eval.content.Avec,...
        'UniformOutput',false);
    As = celltensorconvert(As);
    
    switch s.costfunction
        case {'pathlength coord','acceleration coord'}
            %In the case of these two cost functions, we only care about
            %the connection field, and the metric is always identity.
            %dM is passed a filler value
            [xi,dcost] = arrayfun(@(t,i) get_velocities(t,s,gait,ConnectionEval,As{i},id,1),...
                tpoints,indexList,'UniformOutput',false);
        case {'torque','covariant acceleration','power quality'}
            %In the inertial cases, we need to calculate dM, and the metric
            metrics = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
                samplePoints{:},'spline'),s.metricfield.metric_eval.content.metric,...
                'UniformOutput',false);
            metrics = celltensorconvert(metrics);
            
            dM_set = {};
            for dim = 1:actual_size
                dM_holder = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
                    samplePoints{:},'spline'),s.coriolisfield.coriolis_eval.content.dM{dim},...
                    'UniformOutput',false);
                dM_holder = celltensorconvert(dM_holder);
                dM_set{dim} = dM_holder;
            end
            dMs = {};
            for i = 1:numel(tpoints)
                dMs{i} = {};
                for dim = 1:actual_size
                    dMs{i}{dim} = dM_set{dim}{i};
                end
            end
            
            [xi,dcost] = arrayfun(@(t,i) get_velocities(t,s,gait,ConnectionEval,...
                As{i},metrics{i},dMs{i}),tpoints,indexList,'UniformOutput',false);
        otherwise
            %Otherwise, we're not doing inertial so don't need dM, but we
            %do care about the metric and connection
            metrics = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
                samplePoints{:},'spline'),s.metricfield.metric_eval.content.metric,...
                'UniformOutput',false);
            metrics = celltensorconvert(metrics);
            
            [xi,dcost] = arrayfun(@(t,i) get_velocities(t,s,gait,ConnectionEval,...
                As{i},metrics{i},1),tpoints,indexList,'UniformOutput',false);
    end

	%%%%%%%
	% Integrate cost and displacement into final values
	
	%%
	% Exponential integration for body velocity
	
	% Exponentiate each velocity over the corresponding time step
	expXi = cellfun(@(xi,timestep) se2exp(xi*timestep),xi,num2cell(tsteps),'UniformOutput',false);
	
	% Start off with zero position and displacement
	net_disp_matrix = eye(size(expXi{1}));
	
	% Loop over all the time steps from 1 to n-1, multiplying the
	% transformation into the current displacement
	for i = 1:(length(expXi)-1)
		
		net_disp_matrix = net_disp_matrix * expXi{i};
		
	end
	
	% De-matrixafy the result
	g_theta = atan2(net_disp_matrix(2,1),net_disp_matrix(1,1));
	g_xy = net_disp_matrix(1:2,3);
	
	net_disp_orig = [g_xy;g_theta];
	
	%%
	% Trapezoidal integration for cost
	dcost = [dcost{:}];
	cost = trapz(tpoints,dcost);

end