%Get the resulting fiber motion
function p = se2_integrator_all_terms_fixed(s,p,i,j)

    %get size information about path
    n_dim = size(s.height,3);

    %Prime matrices to hold BVI and displacement data
    p.dZeta = zeros(p.phi_res{i}{j}, n_dim); %discretized body velocity
    p.Zeta = p.dZeta;                                      %BVI
    p.dG = p.dZeta;                                        %discretized global velocity
    p.G = p.dG;                                            %Displacement
	p.dL2 = zeros(p.phi_res{i}{j},1);		   % contribution to pathlength
	p.L = p.dL2;											   % pathlength
	dLoffset = 0;                                          % used in calculating p.dL
	
	p.A_shch_mag = p.dZeta;
	p.A_shch_optimized = p.A_shch_mag;
    

	% Pull out the shape and shape velocity functions
	phi_fun = p.phi_fun{i}{j};
	dphi_fun = p.dphi_fun{i}{j};
	
	% Get the shape and shape derivative at the chosen times
	t = p.time{i}{j};
	t_cell = num2cell(t);
	shape = cellfun(phi_fun,t_cell,'UniformOutput',false);
	dshape = cellfun(dphi_fun,t_cell,'UniformOutput',false);
	
	% Group the individual shape components into cells, so they can be
	% turned into list arguments
	shape_list = cellfun(@(X) num2cell(X),shape,'UniformOutput',false);
	
	
	
	% Evaluate the local connection at each time, in both optimized and
	% original coordinates
	A_interp = @(alpha) cellfun(@(Y) -interpn(s.grid.eval{:},Y,alpha{:}),s.vecfield.eval.content.Avec);
	A_interp_optimized = @(alpha) cellfun(@(Y) -interpn(s.grid.eval{:},Y,alpha{:}),s.vecfield.eval.content.Avec_optimized);

	A_shch = cellfun(A_interp,shape_list,'UniformOutput',false);
	A_shch_optimized = cellfun(A_interp_optimized,shape_list,'UniformOutput',false);
	
	% Get the body velocity at each time step
	%Multiply the connection by the shape diffs
	
	
	A_diff_prod = -sum(A_shch .* repmat(dshape,[n_dim,1]),2);
	A_diff_prod_optimized = sum(A_shch_optimized .* repmat(dshape,[n_dim,1]),2);

	%If there is a ref point connection modifier, calc velocities
	%for it as well
	if isfield(s,'B_ref_point')
		B_diff_prod = sum(s.B_ref_point(p.phi{i,j}(:,2),p.phi{i,j}(:,3)) .* repmat(shape_diff,[n_dim,1]),2);
	end

	%If calculate the pathlength either by the specified metric
	%or by a unit metric (inserted by sys_calcsystem)
	for k = 1:size(shape_diff,1)

		p.dL2(k+dLoffset,j) = shape_diff(k,:) * s.metric(p.phi{i,j}(k,2),p.phi{i,j}(k,3)) * shape_diff(k,:)';

	end
	dLoffset = dLoffset+1;

	%Rearange to get the dZeta values
	for k = 1:n_dim

		p.dZeta(((i-1)*p.phi_res)+1:i*p.phi_res,j,k) = A_diff_prod(1+(k-1)*p.phi_res:k*p.phi_res);
		p.dZeta_optimized(((i-1)*p.phi_res)+1:i*p.phi_res,j,k) = A_diff_prod_optimized(1+(k-1)*p.phi_res:k*p.phi_res);

		p.A_shch_mag(((i-1)*p.phi_res)+1:i*p.phi_res,j,k) = A_shch_mag(1+(k-1)*p.phi_res:k*p.phi_res);
		p.A_shch_optimized_mag(1+(k-1)*p.phi_res:k*p.phi_res,j,k) = A_shch_optimized_mag(1+(k-1)*p.phi_res:k*p.phi_res);

		%If there is a ref point connection modifier, calc velocities
		%for it as well
		if isfield(s,'B_ref_point')
			p.dZeta_refpoint(((i-1)*p.phi_res)+1:i*p.phi_res,j,k) = B_diff_prod(1+(k-1)*p.phi_res:k*p.phi_res);
		end

	end           

    %Use trapezoidal integration
    p.Zeta = cumtrapz(p.dZeta,1);
    p.Zeta_optimized = cumtrapz(p.dZeta_optimized,1);
	if isfield(s,'B_ref_point')
		p.Zeta_refpoint = cumtrapz(p.dZeta_refpoint,1);
	end
	
	%Calculate the pathlength
	p.L = cumtrapz(sqrt(p.dL2),1);
	
	
    %for each path
    for i = 1:n_paths
        
            %get the mean angle during each time step
            mean_theta = [p.Zeta(1,i,3); .5*(p.Zeta(1:end-1,i,3) +p.Zeta(2:end,i,3))];
                
            p.dG(:,i,1) = p.dZeta(:,i,1).*cos(mean_theta)-...
                p.dZeta(:,i,2).*sin(mean_theta);

            p.dG(:,i,2) = p.dZeta(:,i,1).*sin(mean_theta)+...
                p.dZeta(:,i,2).*cos(mean_theta);

            p.dG(:,i,3) = p.dZeta(:,i,3);
            
			%get the mean angle during each time step
            mean_theta_optimized = [p.Zeta_optimized(1,i,3); .5*(p.Zeta_optimized(1:end-1,i,3) +p.Zeta_optimized(2:end,i,3))];
                
            p.dG_optimized(:,i,1) = p.dZeta_optimized(:,i,1).*cos(mean_theta_optimized)-...
                p.dZeta_optimized(:,i,2).*sin(mean_theta_optimized);

            p.dG_optimized(:,i,2) = p.dZeta_optimized(:,i,1).*sin(mean_theta_optimized)+...
                p.dZeta_optimized(:,i,2).*cos(mean_theta_optimized);

            p.dG_optimized(:,i,3) = p.dZeta_optimized(:,i,3);
        
    end
    
%     %take the cumulative sum
%     p.G = cumsum(p.dG,1);

    %Use trapezoidal integration
    p.G = cumtrapz(p.dG,1);
	p.G_optimized = cumtrapz(p.dG_optimized,1);

end