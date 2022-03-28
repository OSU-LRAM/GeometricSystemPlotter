%s = system struct exported from sysplotter
%p = gait struct exported from sysplotter
%costfunction = 'torque' or 'covariant acceleration'
function normalizedPeriod = getNormalizedPeriod(s,gait,costfun)
    
    s.costfunction = costfun;

    t = gait.time_full{1};
    phi_def = gait.phi_def{1}{1}(t);
    t = linspace(0,1,numel(t));
    
    dimension = size(phi_def,2);

    fa=cell(1,dimension);
    % The bounds ensure the fourier series terms have the right period
     % Bound only the frequency to be 2*pi, such that period = 1
    options = fitoptions('fourier4');
    options.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf 2*pi];
    options.Upper = -[-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf -2*pi];
    
    for i=1:1:dimension
        fa{i}=fit(t',phi_def(:,i),'fourier4',options);
    end
    
    nu={'a0';'a1';'b1';'a2';'b2';'a3';'b3';'a4';'b4';'w'};
    
    y = zeros(numel(nu),dimension);
    for i=1:dimension
        for j=1:length(nu)
            y(j,i)=fa{i}.(nu{j});
        end
    end
    
    p = makeGait(y);
    
    [~,cost] = fixed_step_integrator(s,p,[0,1],'interpolated',100);
    
    if strcmpi(s.costfunction,'torque') || strcmpi(s.costfunction,'covariant acceleration') || strcmpi(s.costfunction,'acceleration coord') || strcmpi(s.costfunction,'power quality')
        % With cost as time period, period of the gait is the cost to execute
        % the gait at unit torque squared to the 1/4th power
        normalizedPeriod = cost^(1/4);
    else
        % Drag systems have cost as path length, so no modification is needed
        normalizedPeriod = cost;
    end

end


% function [net_disp_orig, cost] = fixed_step_integrator(s,gait,tspan,ConnectionEval,resolution)
% 
% 	% Duplicate 'resolution' to 'res' if it is a number, or place res at a
% 	% starting resolution if an automatic convergence method is selected
% 	% (automatic convergence not yet enabled)
% 	default_res = 100;
% 	if isnumeric(resolution)
% 		res = resolution;
% 	elseif ischar(resolution) && strcmp(resolution,'autoconverge')
% 		res = default_res;
% 	else
% 		error('Unexpected value for resolution');
% 	end
% 	
% 	% Generate the fixed points from the time span and resolution
% 	tpoints = linspace(tspan(1),tspan(2),res);
% 	tsteps = gradient(tpoints);
%     
%     %Prep interpolation inputs for velocities function
%     shape = zeros(size(s.grid.eval));
%     shape_gait_def_0 = readGait(gait.phi_def,0);
%     actual_size = min(numel(shape),numel(shape_gait_def_0));
%     
%     samplePoints = {};
%     for dim = 1:actual_size
%         samplePoints{dim} = [];
%     end
%     
%     for time = tpoints
%         shape_gait_def = readGait(gait.phi_def,time);
%         shape(1:actual_size) = shape_gait_def(1:actual_size);
%         for dim = 1:numel(shape)
%             samplePoints{dim}(end+1) = shape(dim); 
%         end
%     end
%     
%     indexList = 1:numel(tpoints);
%     id = eye(actual_size);
%     
%     As = cellfun(@(C) -interpn(s.grid.eval{:},C,...
%         samplePoints{:},'spline'),s.vecfield.eval.content.Avec,...
%         'UniformOutput',false);
%     As = celltensorconvert(As);
%     
%     switch s.costfunction
%         case {'pathlength coord','acceleration coord'}
%             %In the case of these two cost functions, we only care about
%             %the connection field, and the metric is always identity.
%             %dM is passed a filler value
%             [xi,dcost] = arrayfun(@(t,i) get_velocities(t,s,gait,ConnectionEval,As{i},id,1),...
%                 tpoints,indexList,'UniformOutput',false);
%         case {'torque','covariant acceleration','power quality'}
%             %In the inertial cases, we need to calculate dM, and the metric
%             metrics = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
%                 samplePoints{:},'spline'),s.metricfield.metric_eval.content.metric,...
%                 'UniformOutput',false);
%             metrics = celltensorconvert(metrics);
%             
%             dM_set = {};
%             for dim = 1:actual_size
%                 dM_holder = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
%                     samplePoints{:},'spline'),s.coriolisfield.coriolis_eval.content.dM{dim},...
%                     'UniformOutput',false);
%                 dM_holder = celltensorconvert(dM_holder);
%                 dM_set{dim} = dM_holder;
%             end
%             dMs = {};
%             for i = 1:numel(tpoints)
%                 dMs{i} = {};
%                 for dim = 1:actual_size
%                     dMs{i}{dim} = dM_set{dim}{i};
%                 end
%             end
%             
%             [xi,dcost] = arrayfun(@(t,i) get_velocities(t,s,gait,ConnectionEval,...
%                 As{i},metrics{i},dMs{i}),tpoints,indexList,'UniformOutput',false);
%         otherwise
%             %Otherwise, we're not doing inertial so don't need dM, but we
%             %do care about the metric and connection
%             metrics = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
%                 samplePoints{:},'spline'),s.metricfield.metric_eval.content.metric,...
%                 'UniformOutput',false);
%             metrics = celltensorconvert(metrics);
%             
%             [xi,dcost] = arrayfun(@(t,i) get_velocities(t,s,gait,ConnectionEval,...
%                 As{i},metrics{i},1),tpoints,indexList,'UniformOutput',false);
%     end
% 
% 	%%%%%%%
% 	% Integrate cost and displacement into final values
% 	
% 	%%
% 	% Exponential integration for body velocity
% 	
% 	% Exponentiate each velocity over the corresponding time step
% 	expXi = cellfun(@(xi,timestep) se2exp(xi*timestep),xi,num2cell(tsteps),'UniformOutput',false);
% 	
% 	% Start off with zero position and displacement
% 	net_disp_matrix = eye(size(expXi{1}));
% 	
% 	% Loop over all the time steps from 1 to n-1, multiplying the
% 	% transformation into the current displacement
% 	for i = 1:(length(expXi)-1)
% 		
% 		net_disp_matrix = net_disp_matrix * expXi{i};
% 		
% 	end
% 	
% 	% De-matrixafy the result
% 	g_theta = atan2(net_disp_matrix(2,1),net_disp_matrix(1,1));
% 	g_xy = net_disp_matrix(1:2,3);
% 	
% 	net_disp_orig = [g_xy;g_theta];
% 	
% 	%%
% 	% Trapezoidal integration for cost
% 	dcost = [dcost{:}];
% 	cost = trapz(tpoints,dcost);
% 
% end

function gait = makeGait(y)

    ndim = size(y,2);
    
    gait.phi_def = cell(1,ndim);
    gait.dphi_def = cell(1,ndim);
    gait.ddphi_def = cell(1,ndim);
    
    for dim = 1:ndim
        
        w = y(end,dim);
        gait.phi_def{dim} = @(t) y(1,dim)+y(2,dim)*cos(w*t)+y(3,dim)*sin(w*t)+y(4,dim)*cos(2*w*t)+...
                                +y(5,dim)*sin(2*w*t)+y(6,dim)*cos(3*w*t)+y(7,dim)*sin(3*w*t)+...
                                +y(8,dim)*cos(4*w*t)+y(9,dim)*sin(4*w*t);
        gait.dphi_def{dim} = @(t) -w*y(2,dim)*sin(w*t)+w*y(3,dim)*cos(w*t)-2*w*y(4,dim)*sin(2*w*t)+...
                                  +2*w*y(5,dim)*cos(2*w*t)-3*w*y(6,dim)*sin(3*w*t)+3*w*y(7,dim)*cos(3*w*t)+...
                                  -4*w*y(8,dim)*sin(4*w*t)+4*w*y(9,dim)*cos(4*w*t);
        gait.ddphi_def{dim} = @(t) -w^2*y(2,dim)*cos(w*t)-w^2*y(3,dim)*sin(w*t)-4*w^2*y(4,dim)*cos(2*w*t)+...
                                   -4*w^2*y(5,dim)*sin(2*w*t)-9*w^2*y(6,dim)*cos(3*w*t)-9*w^2*y(7,dim)*sin(3*w*t)+...
                                   -16*w^2*y(8,dim)*cos(4*w*t)-16*w^2*y(9,dim)*sin(4*w*t);
        
    end
                
end

function state = readGait(gaitFun,t)

    ndim = numel(gaitFun);
    state = zeros(1,ndim);
    
    for i = 1:ndim
        state(i) = gaitFun{i}(t);
    end

end

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

function dcost = torque_cost(M,dM,shape,dshape,ddshape,metric)
% Calculates the incremental cost for an inertial system where cost is torque squared.
% Inputs:
%   M: Mass matrix
%   dM_alphadalpha: Partial of mass matrix with respect to shape variables;
%       must be cell where the ith entry is the partial with respect to the
%       ith shape variable
%   shape: Current shape configuration of system, at which M and
%       dM_alphadalpha were evaluated
%   dshape: Current shape velocity of system
%   ddshape: Current shape acceleration of system

    % Start by calculating the coriolis matrix
    C = calc_coriolis_matrix(dM,shape,dshape);
    % Calculate the torque for this instant of time and return the inner
    % product of the torque with itself
    dtau = M*ddshape(:) + C;
    dcost = dtau'*dtau;
end

function dcost = acceleration_cost(M,dM,shape,dshape,ddshape,metric)
% Calculates the incremental cost for an inertial system where cost is covariant acceleration.
% Inputs:
%   M: Mass matrix
%   dM_alphadalpha: Partial of mass matrix with respect to shape variables;
%       must be cell where the ith entry is the partial with respect to the
%       ith shape variable
%   shape: Current shape configuration of system, at which M and
%       dM_alphadalpha were evaluated
%   dshape: Current shape velocity of system
%   ddshape: Current shape acceleration of system
    C = calc_coriolis_matrix(dM,shape,dshape);
    cov_acc = ddshape(:) + inv(M)*C;
    dcost = cov_acc'*metric*cov_acc;

end

function C = calc_coriolis_matrix(dM,shape,dshape)
    % Start with (dM_alpha/dalpha*qdot)*qdot terms
    C_temp = zeros(length(shape));
    if isa(dM{1},'sym')
        C_temp = sym(C_temp);
    end
    for i = 1:length(dM)
        C_temp = C_temp + dM{i}*dshape(i);
    end
    C = C_temp*dshape(:);
    % Add on the (-1/2)*qdot'*dM_alpha/dalpha*qdot terms
    C_temp = zeros(length(shape),1);
    if isa(dM{1},'sym')
        C_temp = sym(C_temp);
    end
    for i = 1:length(dM)
        C_temp(i) =  -(1/2)*dshape(:)'*dM{i}*dshape(:);
    end
    C = C + C_temp;
    if isa(dM{1},'sym')
        C = simplify(C);
    end
end

function expXi = se2exp(xi)

	% Make sure xi is a column
	xi = xi(:);

	% Special case non-rotating motion
	if xi(3) == 0
		
		expXi = [eye(2) xi(1:2); 0 0 1];
		
	else
		
		z_theta = xi(3);
		
		z_xy = 1/z_theta * [sin(z_theta), 1-cos(z_theta); cos(z_theta)-1, sin(z_theta)] * xi(1:2);
		
		expXi = [ [cos(z_theta), -sin(z_theta); sin(z_theta), cos(z_theta)], z_xy;
			0 0 1];
		
	end


end

function dcost = power_quality_cost(M,dM,shape,dshape,ddshape)
% Calculates the incremental cost for an inertial system where cost is power quality.
% Inputs:
%   M: Mass matrix
%   dM_alphadalpha: Partial of mass matrix with respect to shape variables;
%       must be cell where the ith entry is the partial with respect to the
%       ith shape variable
%   shape: Current shape configuration of system, at which M and
%       dM_alphadalpha were evaluated
%   dshape: Current shape velocity of system
%   ddshape: Current shape acceleration of system

    % Start by calculating the coriolis matrix
    C = calc_coriolis_matrix(dM,shape,dshape);
    % Calculate the torque for this instant of time 
    dtau = M*ddshape(:) + C;
    % Calculate power quality
    dcost = (dshape(:)'*dtau)^2 - ((dshape(:)').^2*dtau.^2);
    dcost = dcost + 100;
end