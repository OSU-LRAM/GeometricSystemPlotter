function p = find_g(s,p)
% get the position change and related data from integrating the shape
% change in p over the system in s

	% Set the initial values of the integrals
	n_dim = size(s.vecfield.eval.content.Avec_optimized,1);
	
	% Iterate over the paths
	n_paths = length(p.phi_fun);
	for i = 1:n_paths
		
		intstart = struct('bvi', {zeros(1,n_dim)},...
					'bvi_opt', {zeros(1,n_dim)},...
					'G', {zeros(1,n_dim)},...
					'G_opt', {zeros(1,n_dim)},...
					'pathlength', {0});

		% Iterate over the segments in that path
		n_segments = length(p.phi_fun{i});
		for j = 1:n_segments
			
			%%%%
			% Integrate the body velocity integral in the original
			% coordinates
			
			
			if ~isfield(p,'fixed_step_integration') || (p.fixed_step_integration == 0)
			
				% Solution for all the displacements
                %odeoptions = odeset('RelTol',1e-6,'AbsTol',1e-12);
                
				sol = ode45(@(t,y) se2_integrator_all_terms...
					(t,y,s,p.phi_fun{i}{j},p.dphi_fun{i}{j}) ...
					,p.time{i}{j}([1 end]) ...
					,zeros(12,1));%,odeoptions);


	% 			% Extract the BVI and displacement functions relative to the start of the segment
	% 			p.bvi_fun_local{i}{j} = @(t) extract_eval(sol,t,1:3);
	% 			p.bvi_opt_fun_local{i}{j} = @(t) extract_eval(sol,t,4:6);
	% 			p.G_fun_local{i}{j} = @(t) extract_eval(sol,t,7:9);
	% 			p.G_opt_fun_local{i}{j} = @(t) extract_eval(sol,t,10:12);

				%%%%%%%%%%%%%%%%%%
				% Evaluate the position and cost loci

				full_solution = deval(sol,p.time{i}{j});
				
			else
			% Use a fixed-step R-K integration
				
				full_solution = se2_integrator_all_terms_fixed(s,p,i,j);
				
				
			end
			
			% With respect to the start of the segment
			p.G_locus{i}{j}.bvi_local = full_solution(1:3,:)';
			p.G_locus{i}{j}.bvi_opt_local = full_solution(4:6,:)';
			p.G_locus{i}{j}.G_local = full_solution(7:9,:)';
			p.G_locus{i}{j}.G_opt_local = full_solution(10:12,:)';
			
			% With respect to the start of the gait
			p.G_locus{i}{j}.bvi = (p.G_locus{i}{j}.bvi_local...
				+repmat(intstart.bvi,size(p.time{i}{j}(:))));
			
			p.G_locus{i}{j}.bvi_opt = (p.G_locus{i}{j}.bvi_opt_local...
				+repmat(intstart.bvi_opt,size(p.time{i}{j}(:))));
			
			x = intstart.G(1);
			y = intstart.G(2);
			theta = intstart.G(3);
			p.G_locus{i}{j}.G = (repmat([x; y; theta],size(p.time{i}{j}(:)')) + ...
								[cos(theta) -sin(theta) 0;
								sin(theta) cos(theta) 0;
								0 0 1]*p.G_locus{i}{j}.G_local')';
			
			x = intstart.G_opt(1);
			y = intstart.G_opt(2);
			theta = intstart.G_opt(3);
			p.G_locus{i}{j}.G_opt = (repmat([x; y; theta],size(p.time{i}{j}(:)')) + ...
								[cos(theta) -sin(theta) 0;
								sin(theta) cos(theta) 0;
								0 0 1] * p.G_locus{i}{j}.G_opt_local')';
							
			% Reset the start values
			intstart.bvi = p.G_locus{i}{j}.bvi(end,:);
			intstart.bvi_opt = p.G_locus{i}{j}.bvi_opt(end,:);
			intstart.G = p.G_locus{i}{j}.G(end,:);
			intstart.G_opt = p.G_locus{i}{j}.G_opt(end,:);

			
			
			%%%%
			% Integrate the pathlength according to the system distance metric
			sol = ode45(@(t,y) pathlength_integrator ...
				(t,y,s,p.phi_fun{i}{j},p.dphi_fun{i}{j}) ...
				,p.time{i}{j}([1 end]) ...
				,zeros(size(intstart.pathlength)));
			
			p.G_locus{i}{j}.S_local =  deval(sol,p.time{i}{j})';
			p.G_locus{i}{j}.S = p.G_locus{i}{j}.S_local+intstart.pathlength;
			
			% reset the starting point
			intstart.pathlength = p.G_locus{i}{j}.S(end);
			
		
		end
		
		%%%%%%%
		% After all segments have been computed, concatenate the segments
		
				
		%%%%%%%%%%%%%%%%%%
		% Evaluate the position and cost loci for the shape change

		p.G_locus_full{i}.bvi = [];
		p.G_locus_full{i}.bvi_opt = [];
		p.G_locus_full{i}.G = [];
		p.G_locus_full{i}.G_opt = [];
		p.G_locus_full{i}.S = [];

		for j = 1:numel(p.G_locus{i})
			p.G_locus_full{i}.bvi = [p.G_locus_full{i}.bvi; p.G_locus{i}{j}.bvi];
			p.G_locus_full{i}.bvi_opt = [p.G_locus_full{i}.bvi_opt; p.G_locus{i}{j}.bvi_opt];
			p.G_locus_full{i}.G = [p.G_locus_full{i}.G; p.G_locus{i}{j}.G];
			p.G_locus_full{i}.G_opt = [p.G_locus_full{i}.G_opt; p.G_locus{i}{j}.G_opt];
			p.G_locus_full{i}.S = [p.G_locus_full{i}.S; p.G_locus{i}{j}.S];
		end
		
		
		
		
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