function Mp = Internal_dissipation_metric_from_curvature_bases(kappa_basis,L,internal_damping_constant)
% Calculate the internal dissipation power metric for a set of curvature bases

	% Specified integration limits
	int_limit = L*[-0.5 0.5];
	
	% Integrate along the body for the internal dissipation metric
	Mp_sol = ode_multistart(@ode45,@(s,Mp) dMetric(s,Mp,kappa_basis,internal_damping_constant),int_limit,int_limit(1),zeros(length(kappa_basis)^2,1));

	Mp = reshape(Mp_sol(int_limit(end)),length(kappa_basis),[]);

end

function dMp = dMetric(s,Mp,kappa_basis,internal_damping_constant) %#ok<INUSL>
% Calculate the local contribution to the internal dissipation power metric

	dMp = zeros(numel(kappa_basis));
	
    for idx1 = 1:size(dMp,1)
        for idx2 = 1:idx1
            
            dMp(idx1,idx2) = kappa_basis{idx1}(s)*kappa_basis{idx2}(s);
            dMp(idx2,idx1) = dMp(idx1,idx2);
            
        end
        
    end
        
    
	dMp = (internal_damping_constant^2) * dMp(:);
	
end