function Mp = LowRE_dissipation_metric_discrete(linklengths,jointangles,L,baseframe,drag_ratio,c)%(geometry,physics,shapeparams)
% Calculate the dissipation power metric for a set of curvature bases

    %Generate backbone geometry and Jacobian from its local definition
    %[h,J,J_full] = backbone(geometry,shapeparams);
    [h, J, J_full] = backbone_from_links(linklengths,jointangles,L,baseframe);

	


% 	% Get the backbone locus, Jacobian, and Local Connection functions
% 	[A, h, J] = LowRE_local_connection(geometry,physics,shapeparams);
    % Get local connection
    [A] = LowRE_local_connection_discrete(linklengths,jointangles,L,baseframe,drag_ratio,c);

    % Get the Metric contribution from each link
    link_metrics = repmat({zeros(size(J_full{1},2))},size(J_full));
    for idx = 1:numel(link_metrics)
        link_metrics{idx} = LowRE_link_metric(J_full{idx},A,h.lengths(idx),drag_ratio,c);
    end

    % Sum the contributions of the links
    Mp = sum(cat(3,link_metrics{:}),3);

    
    
end

function link_metric = LowRE_link_metric(J_full,A,L,drag_ratio,c)
% Calculate the local contribution to the power metric

	% Define the tangential, lateral drag matrix for unit/double drag
	%drag_matrix = [1 0; 0 physics.drag_ratio]*physics.drag_coefficient;
     drag_matrix =     [L      0               0;
                        0    drag_ratio*L       0;
                        0        0           drag_ratio/12*L^3]*c;
                    
                    
    % Calculate the jacobian from shape variables to body velocity of this
    % point on the system:
    
    J_total = J_full * [-A; eye(size(A,2))];
    
    link_metric = J_total' * drag_matrix * J_total;
	
end