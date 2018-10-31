%get shape change loci and, in two dimensions, their projection onto the
%surface-representations of the relevant functions
function p = find_loci(s,p)

    %prime arrays that would otherwise be constructed in loops
	
	% number of paths is the number of shapechange functions specified
    n_paths = numel(p.phi_def);
	
	% number of segments in each path
	n_segments = cellfun(@(x) numel(x),p.phi_def,'UniformOutput',false);

	% number of shape dimensions
    n_dim = numel(s.grid.eval);
    
	% number of position dimensions
	n_g = size(s.dA,1);
	
	% cell array to hold evaluated loci for the paths
    p.phi_locus = cell(size(p.phi_def));
	
	% List the functions on which to project the shape paths
	fields = fieldnames(s);
	project_list = fields(strncmpi('dA',fields,2));	

    %calculate path loci and generate useful forms of the gait function
	for i = 1:n_paths
        
		time_start = 0;
		
		for j = 1:n_segments{i}
			
			
			%%%%%%%%%%%%%%
			% Format time vectors
			
			%generate time vector if necessary
			if length(p.time_def{i}{j}) == 2

				p.time_def{i}{j} = linspace(p.time_def{i}{j}(1), p.time_def{i}{j}(2), p.phi_res{i}{j})';

			else

				% Leave the time_def vector as is 

			end
			
			% Generate the time vector for this segment, starting at zero
			p.time_local{i}{j} = p.time_def{i}{j}-p.time_def{i}{j}(1);
			
			% Generate the time vector for this segment, with zero at the
			% start of the whole path
			p.time{i}{j} = p.time_local{i}{j} + time_start;
			
			%%%%%%%%%%%%%%%
			% shape change function modifications
			
			% shape change definition remapped to start at zero in each
			% segment
			p.phi_fun_local{i}{j} = @(t) p.phi_def{i}{j}(t+p.time_def{i}{j}(1));
			
			% segment change definition relative to overall time in shape
			% change
			p.phi_fun{i}{j} = @(t) p.phi_def{i}{j}(t+p.time_def{i}{j}(1) - time_start);
			
			% shape change velocity definition remapped to start at zero in each
			% segment
			p.dphi_fun_local{i}{j} = @(t) p.dphi_def{i}{j}(t+p.time_def{i}{j}(1));
			
			% segment velocity definition relative to overall time in shape
			% change
			p.dphi_fun{i}{j} = @(t) p.dphi_def{i}{j}(t+p.time_def{i}{j}(1) - time_start);
			
			%%%%%%%%%%%%%%%%
			% evaluate shape change locus			

			%generate path through shape space
			p.phi_locus{i}{j}.shape = p.phi_def{i}{j}(p.time_def{i}{j});
            test1=length(p.phi_locus{i}{j}.shape(1,:));
            test2=length(p.phi_locus{i}{j}.shape(:,1));
            if test1<n_dim
                p.phi_locus{i}{j}.shape=[p.phi_locus{i}{j}.shape,zeros(test2,n_dim-test1)];
            end
            
            if test1>n_dim
                p.phi_locus{i}{j}.shape=p.phi_locus{i}{j}.shape(:,1:n_dim);
            end
            
			%fill in phi_marker
			p.phi_locus{i}{j}.marker.shape = p.phi_marker{i}{j};

			%add z-data for plotting on height functions if applicable
			if n_dim == 2;
				
				% curvature functions
				for k = 1:length(project_list)
					
					for ii = 1:n_g

						p.phi_locus{i}{j}.(project_list{k}){ii,1}...
							= interp2(s.grid.eval{[2 1]},permute(s.(project_list{k}){ii},[2 1])...
							,p.phi_locus{i}{j}.shape(:,1),p.phi_locus{i}{j}.shape(:,2));
						
						if ~isempty(p.phi_marker{i}{j})
							p.phi_locus{i}{j}.marker.(project_list{k}){ii,1}...
								= interp2(s.grid.eval{[2 1]},permute(s.(project_list{k}){ii},[2 1])...
								,p.phi_marker{i}{j}(:,1),p.phi_marker{i}{j}(:,2));
						else
							p.phi_locus{i}{j}.marker.(project_list{k}){ii,1} = [];							
						end
						
					end
					
				end
								
				% coordinate change functions
				for k = 1:n_g
					
					p.phi_locus{i}{j}.Beta{k,1} ...
						= interp2(s.grid.eval{[2 1]},permute(s.B_optimized.eval.Beta{k},[2 1])...
						,p.phi_locus{i}{j}.shape(:,1),p.phi_locus{i}{j}.shape(:,2));

					if ~isempty(p.phi_marker{i}{j})
						
						p.phi_locus{i}{j}.marker.Beta{k,1} ...
							= interp2(s.grid.eval{[2 1]},permute(s.B_optimized.eval.Beta{k},[2 1])...
							,p.phi_marker{i}{j}(:,1),p.phi_marker{i}{j}(:,2));
					else
						p.phi_locus{i}{j}.marker.Beta{k,1} = [];
					end
					
				end
				
			end

			
			% Update the time_start value for the next segment
			time_start = p.time{i}{j}(end);
		end
		
		% Create a full function for this shape change
		start_times = cellfun(@(t) min(t),p.time{i});
		end_times = cellfun(@(t) max(t),p.time{i});
		all_limits = [start_times(:) end_times(:)];
		
		p.phi_fun_full{i} = @(t) concatenate_functions(p.phi_fun{i},all_limits,t);
		
		% Create a full time vector for this shape change 
		p.time_full{i} = cat(1,p.time{i}{:});
		
		% Insert a full locus for this shape change into the locus
		% structure
		p.phi_locus_full{i}.shape = p.phi_fun_full{i}(p.time_full{1});
        
		
		% pull out the appropriate data and concatenate it
		ptemp = cellfun(@(x) x.marker.shape,p.phi_locus{i},'UniformOutput',false);
        test1=length(p.phi_locus_full{i}.shape(1,:));
        test2=length(p.phi_locus_full{i}.shape(:,1));
        if test1<n_dim
            p.phi_locus_full{i}.shape=[p.phi_locus_full{i}.shape,zeros(test2,n_dim-test1)];
        end
        if test1>n_dim
            p.phi_locus_full{i}.shape=p.phi_locus_full{i}.shape(:,1:n_dim);
        end        
		p.phi_locus_full{i}.marker.shape = cat(1,ptemp{:});              

		if n_dim == 2;
			
			% curvature functions
			for k = 1:length(project_list)
				
				ptemp = cellfun(@(x) x.(project_list{k}),p.phi_locus{i},'UniformOutput',false);
				p.phi_locus_full{i}.(project_list{k})...
					= cellfun(@(varargin) cat(1,varargin{:}),ptemp{:},'UniformOutput',false);

				ptemp = cellfun(@(x) x.marker.(project_list{k}),p.phi_locus{i},'UniformOutput',false);
				p.phi_locus_full{i}.marker.(project_list{k})...
					= cellfun(@(varargin) cat(1,varargin{:}),ptemp{:},'UniformOutput',false);

			end

			% coordinate change functions
			ptemp = cellfun(@(x) x.Beta,p.phi_locus{i},'UniformOutput',false);
			p.phi_locus_full{i}.Beta = ...
				cellfun(@(varargin) cat(1,varargin{:}),ptemp{:},'UniformOutput',false);


		end
	end    
end