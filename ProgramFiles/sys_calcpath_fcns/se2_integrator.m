% Function to evaluate displacement and differential cost at each time
function V = se2_integrator(t,X,s,phi_fun,dphi_fun,style,coordinates)

	% Get the shape and shape derivative at the current time
	shape = phi_fun(t);
	shapelist = num2cell(shape);
	dshape = dphi_fun(t);
	
	% Input order to convert ndgrid to meshgrid
	inputorder = [2 1 3:length(shape)];
	
	
	% Get the local connection at the current time, in the new coordinates
	switch coordinates
		
		case 'original'
			
			A = s.A_num(shapelist{:})./s.A_den(shapelist{:});
			
		case  'optimized'
	
			A = zeros(size(s.vecfield.eval.content.Avec_optimized));
			
			for i = 1:size(A,1)
				for j = 1:size(A,2)

					A(i,j) = - interpn(s.grid.eval{:},s.vecfield.eval.content.Avec_optimized{i,j},shapelist{:});

				end
			end
			
		otherwise
			
			error(['Unsupported coordinate definition -' coordinates '- passed to se2_integrator']);
			
	end
	
	% Get the body velocity at the current time
	xi = - A * dshape(:);
		
	% Rotate body velocity into world frame if making a displacement
	% integral
	switch style

		case 'disp'
			
			theta = X(3);
			V = [cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*xi;
			
		case 'bvi'
			
			V = xi;
			
		otherwise
			
			error(['Unsupported style -' style '- passed to se2_integrator']);
			
	end	


end