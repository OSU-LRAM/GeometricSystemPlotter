function [A,Mp,M_full] = LowRE_local_connection_and_dissipation_metric(geometry,physics,shapeparams)
% Calculate the dissipation power metric for a set of curvature bases

    %Generate backbone geometry and its Jacobian from its local definition
    [h,J] = backbone(geometry,shapeparams);

    % Itegrate from one halflength before the midpoint to one halflength after it
    int_limit = geometry.length*[-0.5 0.5];
	
	% Integrate along the body for the power metric
	M_full_sol = ode45(@(s,Mp) LowRE_TotalMetric_infinitesimal...
        (s,h(s),J(s),physics.drag_coefficient,physics.drag_ratio)...
        ,int_limit,zeros(3+numel(shapeparams)));

    % Form the components of the metric into the standard NxN symmetric matrix form
	M_full = reshape(deval(M_full_sol,int_limit(end)),3+numel(shapeparams),[]);
    
    %%%%%%%%%%%
    % Now extract local connection, then use it to pull metric back to
    % shape space
    
    % Pfaffian constraint for Low Reynolds number flow is first three rows
    % of full metric
    pfaffian = M_full(1:3,:);
    
    % Calculate the local connection by multiplying the inverse of the first
    % block in the Pfaffian by the second block
    A = pfaffian(:,1:3)\pfaffian(:,4:end);
    
    % Pull back the full metric onto the shape space
    Mp = [-A' eye(numel(shapeparams))] * M_full * [A; eye(numel(shapeparams))];


end


function dM_full = LowRE_TotalMetric_infinitesimal(s,h,J,c,drag_ratio) %#ok<INUSL>
% Calculate the derivative of the local connection as it's built up along
% the backbone

    % Drag matrix: drag_ratio is the ratio of lateral drag to longitudinal
    % drag, and c is the total drag coefficient.
    drag_matrix = [1 0 0;0 drag_ratio 0;0 0 0]*c;

    % Jacobian from [gcirc rdot]' to motion of point at location s along body
    J_full = [Adjinv(h) TgLginv(h)*J];
    
    % Contribution to the metric from this point is its local drag matrix,
    % pulled back into joint coordinates by its Jacobian.
    dM_full_matrix = J_full.'*drag_matrix*J_full;
    
    % Turn output into a column vector for ode45
    dM_full = dM_full_matrix(:);
    
end