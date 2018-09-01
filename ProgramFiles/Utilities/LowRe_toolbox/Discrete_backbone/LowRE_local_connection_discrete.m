function [A, h, J,Omega] = LowRE_local_connection_discrete(linklengths,jointangles,L,baseframe,drag_ratio,c) %(geometry,physics,shapeparams)
% Calculate the local connection for a set of curvature bases
%
% Inputs:
% geometry: structure defining system geometry
%      geometry.type: how the system geometry is defined 
%         (e.g., links, curvature basis, general curvature)
%      geometry.function: map from shape variables to local backbone
%         deformation (e.g., curvature or joint angles)
%      geometry.length: total length of swimmer
% cparams: value of shape variables

% c: drag per unit length
% drag_ratio: ratio of lateral to longitudinal drag



%Generate backbone geometry and Jacobian from its local definition
%[h,J,J_full] = backbone(geometry,shapeparams);
[h, J, J_full] = backbone_from_links(linklengths,jointangles,L,baseframe);

% Get the Pfaffian contribution from each link
link_force_maps = J_full;
for idx = 1:numel(link_force_maps)
    link_force_maps{idx} = LowRE_body_drag_link(h.pos(idx,:),J_full{idx},h.lengths(idx),drag_ratio,c);
end

% Sum the contributions of the links
omega = sum(cat(3,link_force_maps{:}),3);

% Calculate the local connection by multiplying the inverse of the first
% block in the Pfaffian by the second block
A = omega(:,1:3)\omega(:,4:end);


end


function omega = LowRE_body_drag_link(h,J_full,L,drag_ratio,c)
% Calculate the matrix that maps  

		
	% Local drag on a link.
    %   Longitudinal drag is proportional to longitudinal velocity. 
    %   Lateral drag is proportional to lateral velocity, with a different drag
    %       coefficient, specified by its ratio to the longitudinal drag. 
    %   Rotational drag is proportional to rotational velocity, with a
    %       coefficient based on integrating lateral drag along the length
    %       of the link
    %   c is the absolute drag coefficient
	gcirc_local_to_F_local = ...
        [-L      0               0;
        0    -drag_ratio*L       0;
        0        0           -drag_ratio/12*L^3]*c;
	
    % Transfer force to midpoint-tangent frame by transpose of the
    % adjoint-inverse action
    F_local_to_F_midpoint = transpose(Adjinv(h));
    
    % Linear map from body and shape velocity to force on body frame is
    % product of the two matrices calculated above and the full Jacobian
    % for the link
    omega = F_local_to_F_midpoint * gcirc_local_to_F_local * J_full;
	

end