function [A,Mp,M_full] = LowRE_connection_and_metric_discrete(linklengths,shapeparams,L,baseframe,drag_ratio,c)%(geometry,physics,shapeparams)
% Calculate the dissipation power metric for a set of curvature bases

    %Generate backbone geometry and Jacobian from its local definition
    %[h,J,J_full] = backbone(geometry,shapeparams);
    [h, J, J_full] = backbone_from_links(linklengths,shapeparams,L,baseframe);

	

% 	% Get the backbone locus, Jacobian, and Local Connection functions
% 	[A, h, J] = LowRE_local_connection(geometry,physics,shapeparams);
    % Get local connection
    [A] = LowRE_local_connection_discrete(linklengths,shapeparams,L,baseframe,drag_ratio,c);

    % Get the Metric contribution from each link
    link_metrics = repmat({zeros(size(J_full{1},2))},size(J_full));
    for idx = 1:numel(link_metrics)
        link_metrics{idx} = LowRE_link_metric_full(J_full{idx},A,h.lengths(idx),drag_ratio,c);
    end

    % Sum the contributions of the links
    M_full = sum(cat(3,link_metrics{:}),3);

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

function link_metric_full = LowRE_link_metric_full(J_full,A,L,drag_ratio,c)
% Calculate the local contribution to the power metric


	% Define the tangential, lateral drag matrix for unit/double drag
	%drag_matrix = [1 0; 0 physics.drag_ratio]*physics.drag_coefficient;
     drag_matrix =     [L      0               0;
                        0    drag_ratio*L       0;
                        0        0           drag_ratio/12*L^3]*c;

    % Calculate the jacobian from shape variables to body velocity of this
    % point on the system:
        
    link_metric_full = J_full' * drag_matrix * J_full;
	
end