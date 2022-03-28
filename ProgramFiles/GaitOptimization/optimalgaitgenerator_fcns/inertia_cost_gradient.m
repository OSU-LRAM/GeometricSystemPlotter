function cost_grad = inertia_cost_gradient(s,n,y,g,gait,EvaluationMethod)
% Calculates the gradient of cost for inertial systems.
% Inputs:
%   s: System struct used by sysplotter.
%   n: Number of points at which gait should be evaluated in the shape
%       space.
%   y: Fourier coefficients that parametrize the gait.
%   g: Time period over which gait is executed.
%   gait: Struct containing fields:
%       phi_def: array function that returns shape at time value t
%       dphi_def: array function that returns shape velocity at time t
%       ddphi_def: array function that returns shape acceleration at time t
%   EvaluationMethod: String representing how the gradient of cost should
%   be calculated; provide 'discrete' to evaluate at 100 discrete values or
%   'ode45' if you would like the gradient of cost to be integrated using
%   ode45. Other values will result in error.

    % Contribution to gradient from the movement of each point due to
    % change in fourier coefficients
    [grad_alphaddot,grad_alphadot,grad_alpha] = shape_grad(n,y,g);

    cost_grad = zeros(size(grad_alpha));
    tspan = [0 g];
    if strcmpi(EvaluationMethod,'discrete')
        num_pts = 100;
        t_pts = linspace(0,g,num_pts);
        del_t = t_pts(2) - t_pts(1);
        
        %Prep interpolation inputs for gradient calcs
        shape = zeros(size(s.grid.eval));
        shape_gait_def_0 = readGait(gait.phi_def,0);
        actual_size = min(numel(shape),numel(shape_gait_def_0));
        
        samplePoints = {};
        for dim = 1:actual_size
            samplePoints{dim} = [];
        end
        
        for time = t_pts
            shape_gait_def = readGait(gait.phi_def,time);
            shape(1:actual_size) = shape_gait_def(1:actual_size);
            for dim = 1:numel(shape)
                samplePoints{dim}(end+1) = shape(dim);
            end
        end
        
        %Batch interpolate the metric at each point along the gait
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
        for i = 1:numel(t_pts)
            dMs{i} = {};
            for dim = 1:actual_size
                dMs{i}{dim} = dM_set{dim}{i};
            end
        end
        
        empty_ddM = cell(size(s.coriolisfield.coriolis_eval.content.ddM));
        ddM_set = empty_ddM;
        for dim = 1:numel(ddM_set)
            ddM_holder = cellfun(@(C) interpn(s.grid.metric_eval{:},C,...
                samplePoints{:},'spline'),s.coriolisfield.coriolis_eval.content.ddM{dim},...
                'UniformOutput',false);
            ddM_holder = celltensorconvert(ddM_holder);
            ddM_set{dim} = ddM_holder;
        end
        ddMs = {};
        for i = 1:numel(t_pts)
            ddMs{i} = empty_ddM;
            for dim = 1:numel(empty_ddM)
                ddMs{i}{dim} = ddM_set{dim}{i};
            end
        end
                
        switch s.costfunction
            case 'torque'
                for k = 1:length(t_pts)
                    del_cost = inertia_gradient_helper(t_pts(k),[],s,gait,grad_alpha,grad_alphadot,grad_alphaddot,metrics{k},dMs{k},ddMs{k});
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
            case 'covariant acceleration'
                for k = 1:length(t_pts)
                    del_cost = acceleration_gradient_helper(t_pts(k),[],s,gait,grad_alpha,grad_alphadot,grad_alphaddot,metrics{k},dMs{k},ddMs{k});
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
            case 'acceleration coord'
                for k = 1:length(t_pts)
                    del_cost = accelerationcoord_gradient_helper(t_pts(k),[],s,gait,grad_alphaddot,metrics{k},dMs{k},ddMs{k});
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
            case 'power quality'
                for k = 1:length(t_pts)
                    del_cost = powerquality_gradient_helper(t_pts(k),[],s,gait,grad_alpha,grad_alphadot,grad_alphaddot,metrics{k},dMs{k},ddMs{k});
                    cost_grad = cost_grad + reshape(del_cost,size(cost_grad)).*del_t;
                end
        end
        % Reset gradient of fourier frequency to be zero to prevent changes
        % to it
        cost_grad(end,:) = 0;
    elseif strcmpi(EvaluationMethod,'ode45')
        switch s.costfunction
            case 'torque'
                sol = ode45(@(t,y) inertia_gradient_helper(t,y,s,gait,grad_alpha,grad_alphadot,grad_alphaddot),tspan,cost_grad);
            case 'covariant acceleration'
                sol = ode45(@(t,y) acceleration_gradient_helper(t,y,s,gait,grad_alpha,grad_alphadot,grad_alphaddot),tspan,cost_grad);
            case 'acceleration coord'
                sol = ode45(@(t,y) accelerationcoord_gradient_helper(t,y,s,gait,grad_alphaddot),tspan,cost_grad);
            case 'power quality'
                sol = ode45(@(t,y) powerquality_gradient_helper(t,y,s,gait,grad_alpha,grad_alphadot,grad_alphaddot),tspan,cost_grad);
        end
        % Extract the final motion
        cost_grad = reshape(deval(sol,tspan(end)),size(cost_grad));
        % Reset gradient of fourier frequency to be zero to prevent changes
        % to it
        cost_grad(end,:) = 0;
    else
        error('Untenable option provided for EvaluationMethod!')
    end
end