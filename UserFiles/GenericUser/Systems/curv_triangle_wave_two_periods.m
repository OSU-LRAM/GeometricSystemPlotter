% (516189106917469267*sign((a1*cos(4*s*pi) + a2*sin(4*s*pi))/(a1^2 + a2^2)^(1/2))*abs(a1*cos(4*s*pi) + a2*sin(4*s*pi))^20*(a1^2 + a2^2)^(1/2))/(144115188075855872*abs(a1^2 + a2^2)^10)
function output = curv_triangle_wave_two_periods(params,mode)

% Turn params into a cell matrix
params = num2cell(params);

switch mode
    
    case 'curvature'

		output = @(s) curv_fun(s,params{:});


    case 'orientation'
        
		%% Padded length of unit backbone
		all_limits = [-.51 0 .51];
		
		%% Make dummy integration function
		curv_fun_dummy = curv_triangle_wave_two_periods(cell2mat(params),'curvature');
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
		d_curv_dp_fun_dummy = curv_triangle_wave_two_periods(cell2mat(params),'dcurvature');
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
	t2 = s.*pi.*4.0;
	t3 = a1.^2;
	t4 = a2.^2;
	t5 = t3+t4;
	t6 = cos(t2);
	t7 = a1.*t6;
	t8 = sin(t2);
	t9 = a2.*t8;
	t10 = t7+t9;
	t11 = abs(t10);
	t12 = t11.^2;
	t13 = t12.^2;
	t14 = t13.^2;
	t15 = t14.^2;
	out1 = sqrt(t5).*t13.*t15.*sign(1.0./sqrt(t5).*t10).*1.0./abs(t5).^10.*3.581781447252944;
end

