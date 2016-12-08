function output = shchf_butterfly(input_mode,pathnames)
 %[r1,r2] = convert.old_to_new_points(alpha1,alpha2);
	% Default argument
	if ~exist('input_mode','var')
		
		input_mode = 'initialize';
		
	end
		
	% Name the .mat file with the fourier coefficients
	fourierfile = 'butterfly_params';
	
	
	switch input_mode
		
		case 'name'
			
			output = 'Butterfly gait';
			
		case 'dependency'
			
			output.dependency = {};
			
		case 'initialize'

			%%%%
			%%
			%Path definitions

			%Fourier Coefficients
			% A is cosine, B is sine, number is which joint angle
			load(fourierfile)
			

			p.phi_def{1} = {@(t) [fourier_eval(t(:)+pi,A1,B1,2*pi);fourier_eval(t(:)+pi,A2,B2,2*pi)]'};
			
			
			%marker locations
			p.phi_marker = [];
			
			%arrows to plot
			p.phi_arrows = {{2}};

			%time to run path
			p.time_def{1} = {[0 2*pi]};


			%path resolution
			p.phi_res{1} = {50};
			



			%%%%
			%Save the shch properties
			output = p;
	end
	
end