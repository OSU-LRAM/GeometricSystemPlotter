function [cBVI, cBVI_opt] = integrate_cBVI(s,p)
% Integrate up the corrected body velocity integral

	% number of dimensions for the cBVI
	n_g = size(s.DA,1);
	
	cBVI = cell(size(p.phi_fun_full));
	cBVI_opt = cBVI;

	for i = 1:numel(p.phi_fun_full)
		
		% Extract the evaluated shape
		X = p.phi_locus_full{i}.shape;
		
		% remove any duplicate values
		[junk,ix] = unique(X,'rows');
		X = X(sort(ix),:);
		
		
		% Simple integration is just to throw this at doubleintegral and
		% polyorient. Try this first, unless user has specified otherwise
		if strcmp(p.cBVI_method{i}{1},'simple')
			
			try
				
				cBVI{i} = polysign(X(:,1),X(:,2)) * arrayfun(@(j)...
					doubleintegral(...
					@(a,b) interp2(s.grid.eval{2},s.grid.eval{1},...
					permute(s.DA{j}, [2 1]),a,b,'cubic')...
					,struct('type','polygon','x',X(:,1),'y',X(:,2))...
					,struct('method','gauss','tol',1e-6))...
					,(1:n_g)');

				cBVI_opt{i} = polysign(X(:,1),X(:,2)) * arrayfun(@(j)...
					doubleintegral(...
					@(a,b) interp2(s.grid.eval{2},s.grid.eval{1},...
					permute(s.DA_optimized{j}, [2 1]),a,b,'cubic')...
					,struct('type','polygon','x',X(:,1),'y',X(:,2))...
					,struct('method','gauss','tol',1e-6))...
					,(1:n_g)');
				
            catch ME
				
                if numel(s.grid.eval) == 2
                    warning('Gait curve not sufficiently convex for simple integration. Specify either ''delaunay'' or ''section'' in p.complex_curve{i} to invoke that integration method')
                end
			end
			
		else
			
			switch p.cBVI_method{i}{1}
				
				case 'delaunay'
					
					% make a delaunay triangulation and integrate over it.
					% note that this slows down rapidly with the number of
					% vertices

					% Polygon edges, for discarding triangles outside of
					% the gait
					C = [(1:size(X,1))' [2:size(X,1) 1]' ];

					% Delaunay triangulation and a list of triangles inside
					% the polygon
					warning off
					dt = DelaunayTri(X(:,1), X(:,2), C);
					warning on
					io = dt.inOutStatus();
					
					% Extract out the triangulation, and update the polygon
					% with any points merged or generated in the
					% triangulation
					Triangulation = dt.Triangulation(io,:);
					X = dt.X;
					
					% Get the signed area integral over each triangle, and
					% sum them
					cBVI{i} = arrayfun(@(j) ...
						sum(...
						tri_orient_check(Triangulation) ...
						.* arrayfun(@(k)...
						doubleintegral(...
						@(a,b) interp2(s.grid.eval{2},s.grid.eval{1},...
						permute(s.DA{j}, [2 1]),a,b,'cubic')...
						,struct('type','polygon','x',X(Triangulation(k,:),1)...
						,'y',X(Triangulation(k,:),2))...
						,struct('method','gauss','tol',1e-6)),(1:length(Triangulation))'),1)...
						,(1:n_g)');

					cBVI_opt{i} = arrayfun(@(j) ...
						sum(...
						tri_orient_check(Triangulation) ...
						.* arrayfun(@(k)...
						doubleintegral(...
						@(a,b) interp2(s.grid.eval{2},s.grid.eval{1},...
						permute(s.DA_optimized{j}, [2 1]),a,b,'cubic')...
						,struct('type','polygon','x',X(Triangulation(k,:),1)...
						,'y',X(Triangulation(k,:),2))...
						,struct('method','gauss','tol',1e-6)),(1:length(Triangulation))'),1)...
						,(1:n_g)');
					
				case 'section'
					
					% Apply the simple integration rule over each section
					% of the gait, then sum them
					cBVI_section = zeros(n_g,numel(p.phi_fun{i}));
					for j = 1:numel(p.phi_fun{i})
					
						% Extract the evaluated shape
						X = p.phi_locus{i}{j}.shape;
						% remove any duplicate values
						[junk,ix] = unique(X,'rows');
						X = X(sort(ix),:);

						

						try

							cBVI_section_original(:,j) = polysign(X(:,1),X(:,2)) * arrayfun(@(k)...
								doubleintegral(...
								@(a,b) interp2(s.grid.eval{2},s.grid.eval{1},...
								permute(s.DA{k}, [2 1]),a,b,'cubic')...
								,struct('type','polygon','x',X(:,1),'y',X(:,2))...
								,struct('method','gauss','tol',1e-6))...
								,(1:n_g)');
							
							cBVI_section_optimized(:,j) = polysign(X(:,1),X(:,2)) * arrayfun(@(k)...
								doubleintegral(...
								@(a,b) interp2(s.grid.eval{2},s.grid.eval{1},...
								permute(s.DA_optimized{k}, [2 1]),a,b,'cubic')...
								,struct('type','polygon','x',X(:,1),'y',X(:,2))...
								,struct('method','gauss','tol',1e-6))...
								,(1:n_g)');

						catch

							warning('Gait section not sufficiently convex for simple integration. Specify ''delaunay'' integration for this gait')

						end
						
					end
					
					cBVI{i} = sum(cBVI_section_original,2);
					cBVI_opt{i} = sum(cBVI_section_optimized,2);

			
			end
			
		end
		
		
		
		
	end

end