% (1032378213834938331*sign((a1*cos(2*s*pi) + a2*sin(2*s*pi))/(a1^2 + a2^2)^(1/2))*abs(a1*cos(2*s*pi) + a2*sin(2*s*pi))^20*(a1^2 + a2^2)^(1/2))/(288230376151711744*abs(a1^2 + a2^2)^10)
function output = curv_triangle_wave_one_period(params,mode)

% Turn params into a cell matrix
params = num2cell(params);

switch mode
    
    case 'curvature'

		output = @(s) curv_fun(s,params{:});


    case 'orientation'
        
		%% Padded length of unit backbone
		all_limits = [-.51 0 .51];
		
		%% Make dummy integration function
		curv_fun_dummy = curv_triangle_wave_one_period(cell2mat(params),'curvature');
		curvature = @(s,~) curv_fun_dummy(s);
		
		%% Integral of the integrand function along s
		output = ode_multistart(@ode45,curvature,all_limits,0,0);

        
    case 'dcurvature'
        
		%% Create a dummy function that takes in a vector of parameters
		%% including s, and deals them into individual function parameters
		curv_intermediate = @(all_params) vector_to_list_input(@curv_fun,all_params);
		
		%% Create a function that takes the jacobian of the dummy function
		fulljacobian = @(s) jacobianest(curv_intermediate,[s,params{:}]);
		
		%% Create a function that truncates the s-derivative from the full jacobian
		output = @(s) reshape_truncate_jacobian(fulljacobian(s));

            
    case 'dcurvature_int'
        
		%% Padded length of unit backbone
		all_limits = [-.51 0 .51];
		
		%% Make dummy integration function
		d_curv_dp_fun_dummy = curv_triangle_wave_one_period(cell2mat(params),'dcurvature');
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

function out1 = curv_fun(s,a1,a2)
	t2 = a1.^2;
	t3 = a2.^2;
	t4 = s.*pi.*2.0;
	t5 = cos(t4);
	t6 = sin(t4);
	t9 = t2+t3;
	t7 = a1.*t5;
	t8 = a2.*t6;
	t10 = t7+t8;
	out1 = sqrt(t9).*sign(1.0./sqrt(t9).*t10).*1.0./abs(t9).^10.*abs(t10).^20.*3.581781447252943;
	end
end

