%This function makes it so curvature definition files can be as simple as
%possible.  It handles all the different ways each curve definition could
%possible be called in our functionality
function output = curvEvalModeHandler(params,mode,curve_fun)

% Turn params into a cell matrix
params = num2cell(params);

%How do we want to use this curvature function
switch mode
    
    %Just return the function handle to the curvature function
    case 'curvature'

		output = @(s) curve_fun(s,params{:});


    %Get orientation of the salp along the backbone givent the curvature
    %function
    case 'orientation'
        
		%% Padded length of unit backbone
		all_limits = [-.51 0 .51];
		
		%% Make dummy integration function
		curv_fun_dummy = curvEvalModeHandler(cell2mat(params),'curvature',curve_fun);
		curvature = @(s,~) curv_fun_dummy(s);
		
		%% Integral of the integrand function along s
		output = ode_multistart(@ode45,curvature,all_limits,0,0);

        
    %Get curvature rate of change for doing derivative shenanigans
    case 'dcurvature'
        
		%% Create a dummy function that takes in a vector of parameters
		%% including s, and deals them into individual function parameters
		curv_intermediate = @(all_params) vector_to_list_input(curve_fun,all_params);
		
		%% Create a function that takes the jacobian of the dummy function
		fulljacobian = @(s) jacobianest(curv_intermediate,[s,params{:}]);
		
		%% Create a function that truncates the s-derivative from the full jacobian
		output = @(s) reshape_truncate_jacobian(fulljacobian(s));

    %Initialization function for jacobian shenanigans
    case 'dcurvature_int'
        
		%% Padded length of unit backbone
		all_limits = [-.51 0 .51];
		
		%% Make dummy integration function
		d_curv_dp_fun_dummy = curvEvalModeHandler(cell2mat(params),'dcurvature',curve_fun);
		dcurvature = @(s,~) d_curv_dp_fun_dummy(s)';
		
		%% Integral of the integrand function along s
		dummy_output = ode_multistart(@ode45,dcurvature,all_limits,0,zeros(size(params(:).')));
		output = @(t) transpose(dummy_output(t));

                                   
end

end

function output = vector_to_list_input(funhandle,all_params)

    all_params = num2cell(all_params);
    
    output = funhandle(all_params{:});
    
end


function output = reshape_truncate_jacobian(J)

    output = J(:,2:end);
    
end