function p = find_g(s,p)
% get the position change and related data from integrating the shape
% change in p over the system in s

	% Set the initial values of the integrals
	n_dim = size(s.height,1);
	
	% Iterate over the paths
	n_paths = length(p.phi_fun);
	for i = 1:n_paths
		
		intstart = struct('bvi', {zeros(n_dim,1)},...
					'bvi_opt', {zeros(n_dim,1)},...
					'G', {zeros(n_dim,1)},...
					'G_opt', {zeros(n_dim,1)},...
					'pathlength', {0});

		% Iterate over the segments in that path
		n_segments = length(p.phi_fun{i});
		for j = 1:n_segments
			
			%%%%
			% Integrate the body velocity integral in the original
			% coordinates
			
			% Solution for all the displacements
			sol = ode45(@(t,y) se2_integrator_all_terms...
 				(t,y,s,p.phi_fun{i}{j},p.dphi_fun{i}{j}) ...
				,p.time{i}{j}([1 end]) ...
				,zeros(12,1));
			

			% Extract the BVI and displacement functions relative to the start of the segment
			p.bvi_fun_local{i}{j} = @(t) extract_eval(sol,t,1:3);
			p.bvi_opt_fun_local{i}{j} = @(t) extract_eval(sol,t,4:6);
			p.G_fun_local{i}{j} = @(t) extract_eval(sol,t,7:9);
			p.G_opt_fun_local{i}{j} = @(t) extract_eval(sol,t,10:12);
			
			% Ditto for the functions relative to the start of the shape
			% change
			
			% bvi original
			p.bvi_fun{i}{j} = @(t) (repmat(intstart.bvi,size(t(:)')) + p.bvi_fun_local{i}{j}(t)')';
			intstart.bvi = p.bvi_fun{i}{j}(p.time{i}{j}(end))';
			
			% bvi optimized
			p.bvi_opt_fun{i}{j} = @(t) (repmat(intstart.bvi_opt,size(t(:)'))...
				+ p.bvi_opt_fun_local{i}{j}(t)')';
			intstart.bvi_opt = p.bvi_opt_fun{i}{j}(p.time{i}{j}(end))';
			
			% displacement original
			x = intstart.G(1);
			y = intstart.G(2);
			theta = intstart.G(3);
			p.G_fun{i}{j} = @(t) (repmat([x; y; theta],size(t(:)')) + ...
								[cos(theta) -sin(theta) 0;
								sin(theta) cos(theta) 0;
								0 0 1] * p.G_fun_local{i}{j}(t)')';
			intstart.G = p.G_fun{i}{j}(p.time{i}{j}(end))';
			
			% displacement optimized
			x = intstart.G_opt(1);
			y = intstart.G_opt(2);
			theta = intstart.G_opt(3);
			p.G_opt_fun{i}{j} = @(t) (repmat([x; y; theta],size(t(:)')) + ...
								[cos(theta) -sin(theta) 0;
								sin(theta) cos(theta) 0;
								0 0 1] * p.G_opt_fun_local{i}{j}(t)')';

			% set the initial value for the next segment
			intstart.G_opt = p.G_opt_fun{i}{j}(p.time{i}{j}(end))';			
			
			
			
% 			sol = ode45(@(t,y) se2_integrator...
% 				(t,y,s,p.phi_fun{i}{j},p.dphi_fun{i}{j},'bvi','original') ...
% 				,p.time{i}{j}([1 end]) ...
% 				,zeros(size(intstart.bvi)));
% 			
% 			% Extract the BVI function relative to the start of the segment
% 			p.bvi_fun_local{i}{j} = @(t) deval(sol,t)';
% 			
% 			% Extract the BVI function relative to the start of the shape
% 			% change
% 			p.bvi_fun{i}{j} = @(t) (repmat(intstart.bvi,size(t(:)')) + p.bvi_fun_local{i}{j}(t)')';
% 			
% 			% set the initial value for the next segment
% 			intstart.bvi = p.bvi_fun{i}{j}(p.time{i}{j}(end))';
% 			
% 			%%%%
% 			% Integrate the body velocity integral in the optimized
% 			% coordinates
% 			sol = ode45(@(t,y) se2_integrator...
% 				(t,y,s,p.phi_fun{i}{j},p.dphi_fun{i}{j},'bvi','optimized') ...
% 				,p.time{i}{j}([1 end]) ...
% 				,zeros(size(intstart.bvi_opt)));
% 			
% 			% Extract the BVI function relative to the start of the segment
% 			p.bvi_opt_fun_local{i}{j} = @(t) deval(sol,t)';
% 			
% 			% Extract the BVI function relative to the start of the shape
% 			% change
% 			p.bvi_opt_fun{i}{j} = @(t) (repmat(intstart.bvi_opt,size(t(:)'))...
% 				+ p.bvi_opt_fun_local{i}{j}(t)')';
% 
% 			% set the initial value for the next segment
% 			intstart.bvi_opt = p.bvi_opt_fun{i}{j}(p.time{i}{j}(end))';
% 			
% 			%%%%
% 			% Integrate the position trajectory in the original
% 			% coordinates
% 			sol = ode45(@(t,y) se2_integrator...
% 				(t,y,s,p.phi_fun{i}{j},p.dphi_fun{i}{j},'disp','original') ...
% 				,p.time{i}{j}([1 end]) ...
% 				,zeros(size(intstart.G)));
% 			
% 			% Extract the position trajectory wrt the start of this segment
% 			p.G_fun_local{i}{j} = @(t) deval(sol,t)';
% 
% 			% Extract the position trajectory wrt the start of the shape
% 			% change
% 			x = intstart.G(1);
% 			y = intstart.G(2);
% 			theta = intstart.G(3);
% 			p.G_fun{i}{j} = @(t) (repmat([x; y; theta],size(t(:)')) + ...
% 								[cos(theta) -sin(theta) 0;
% 								sin(theta) cos(theta) 0;
% 								0 0 1] * p.G_fun_local{i}{j}(t)')';
% 			
% 			% set the initial value for the next segment
% 			intstart.G = p.G_fun{i}{j}(p.time{i}{j}(end))';
% 			
% 			%%%%
% 			% Integrate the position trajectory in the optimal coordinates
% 			sol = ode45(@(t,y) se2_integrator...
% 				(t,y,s,p.phi_fun{i}{j},p.dphi_fun{i}{j},'disp','optimized') ...
% 				,p.time{i}{j}([1 end]) ...
% 				,zeros(size(intstart.G_opt)));
% 			
% 			% Extract the position trajectory wrt the start of this segment
% 			p.G_opt_fun_local{i}{j} = @(t) deval(sol,t)';
% 
% 			% Extract the position trajectory wrt the start of the shape
% 			% change
% 			x = intstart.G_opt(1);
% 			y = intstart.G_opt(2);
% 			theta = intstart.G_opt(3);
% 			p.G_opt_fun{i}{j} = @(t) (repmat([x; y; theta],size(t(:)')) + ...
% 								[cos(theta) -sin(theta) 0;
% 								sin(theta) cos(theta) 0;
% 								0 0 1] * p.G_opt_fun_local{i}{j}(t)')';
% 
% 			% set the initial value for the next segment
% 			intstart.G_opt = p.G_opt_fun{i}{j}(p.time{i}{j}(end))';			
			
			%%%%
			% Integrate the pathlength according to the system distance metric
			sol = ode45(@(t,y) pathlength_integrator ...
				(t,y,s.metric,p.phi_fun{i}{j},p.dphi_fun{i}{j}) ...
				,p.time{i}{j}([1 end]) ...
				,zeros(size(intstart.pathlength)));

			% Extract the cost of this segment
			p.S_fun_local{i}{j} = @(t) deval(sol,t)';
			
			% Extract the cost with respect to the start of the shape
			% change
			p.S_fun{i}{j} = @(t) (intstart.pathlength + p.S_fun_local{i}{j}(t)')';
			

			% set the initial value for the next segment
			intstart.pathlength = p.S_fun{i}{j}(p.time{i}{j}(end))';

			%%%%%%%%%%%%%%%%%%
			% Evaluate the position and cost loci
			p.G_locus{i}{j}.bvi = p.bvi_fun{i}{j}(p.time{i}{j});
			p.G_locus{i}{j}.bvi_opt = p.bvi_opt_fun{i}{j}(p.time{i}{j});
			p.G_locus{i}{j}.G = p.G_fun{i}{j}(p.time{i}{j});
			p.G_locus{i}{j}.G_opt = p.G_opt_fun{i}{j}(p.time{i}{j});
			p.G_locus{i}{j}.S = p.S_fun{i}{j}(p.time{i}{j});
			
			

		end
		
		%%%%%%%
		% After all segments have been computed, concatenate the segments
		
		% Get the time limits
		start_times = cellfun(@(t) min(t),p.time{i});
		end_times = cellfun(@(t) max(t),p.time{i});
		all_limits = [start_times(:) end_times(:)];
		
		p.bvi_fun_full{i} = @(t) concatenate_functions(p.bvi_fun{i},all_limits,t);
		p.bvi_opt_fun_full{i} = @(t) concatenate_functions(p.bvi_opt_fun{i},all_limits,t);

		p.G_fun_full{i} = @(t) concatenate_functions(p.G_fun{i},all_limits,t);
		p.G_opt_fun_full{i} = @(t) concatenate_functions(p.G_opt_fun{i},all_limits,t);

		p.S_fun_full{i} = @(t) concatenate_functions(p.S_fun{i},all_limits,t);
		
		%%%%%%%%%%%%%%%%%%
		% Evaluate the position and cost loci for the shape change
		p.G_locus_full{i}.bvi = p.bvi_fun_full{i}(p.time_full{i});
		p.G_locus_full{i}.bvi_opt = p.bvi_opt_fun_full{i}(p.time_full{i});
		p.G_locus_full{i}.G = p.G_fun_full{i}(p.time_full{i});
		p.G_locus_full{i}.G_opt = p.G_opt_fun_full{i}(p.time_full{i});
		p.G_locus_full{i}.S = p.S_fun_full{i}(p.time_full{i});
		
		
		
		
		
	end

	% Get the cBVI for the system
	if ~strcmp(p.cBVI_method{i}{1},'none')
		[p.cBVI p.cBVI_opt] = integrate_cBVI(s,p);
	end

end

function v = extract_eval(sol,t,vals)

	 u = deval(sol,t)';
	 v = u(:,vals);

end