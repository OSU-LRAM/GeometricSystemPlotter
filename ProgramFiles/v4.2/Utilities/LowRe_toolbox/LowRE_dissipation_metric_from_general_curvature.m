function Mp = LowRE_dissipation_metric_from_general_curvature(curvdef,cparams,L,c,drag_ratio)
% Calculate the dissipation power metric for a set of curvature bases

	% Specified integration limits
	int_limit = L*[-0.5 0.5];
	
	% Define the tangential, lateral drag matrix for unit/double drag
	drag = [1 0; 0 drag_ratio]*c;

	% Get the backbone locus, Jacobian, and Local Connection functions
	[A, h, J,Omega] = LowRE_local_connection_from_general_curvature(curvdef,cparams,L,c,drag_ratio);

	% Integrate along the body for the power metric
	Mp_sol = ode_multistart(@ode45,@(s,Mp) dMetric(s,Mp,A,h,J,drag),int_limit,int_limit(1),zeros(length(cparams)^2,1));

	Mp = reshape(Mp_sol(int_limit(end)),length(cparams),[]);

end


function localJ = local_body_velocity_J(A,h,J)
% Calculate the Jacobian from shape parameter velocity to local tangential
% and normal velocity

	R = [cos(h(3)) sin(h(3));
		-sin(h(3)) cos(h(3))];

	%Velocity is sum of body velocity and velocity within body frame
	localJ = R*(-A(1:2,:) + J(1:2,:) + [-h(2)*(-A(3,:)); h(1)*(-A(3,:))]);

end

function dMp = dMetric(s,Mp,A,h,J,drag) %#ok<INUSL>
% Calculate the local contribution to the power metric

	localJ = local_body_velocity_J(A,h(s),J(s));
	
	dMp = localJ.'*drag*localJ;
	
	dMp = dMp(:);
	
end