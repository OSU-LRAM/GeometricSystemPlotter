function Mp = Inertial_metric_discrete(geometry,physics,jointangles)
% Calculate the power dissipation metric for an n-link chain under a
% viscous resitive-force model
%
% Inputs:
%
%   geometry: Structure containing information about geometry of chain.
%       Fields required for this function are:
%
%       linklengths: A vector of link lengths defining a kinematic chain
%
%       baseframe (optional): Specification of which link on the body should be
%           treated as the base frame. Default behavior is to put the base
%           frame at the center of the chain, either at the midpoint of the
%           middle link, or at the joint between the two middle links and at
%           the mean of their orientation. Options for this field are
%           specified in the documentation of N_link_chain.
%
%       length (optional): Total length of the chain. If specified, the elements of
%           will be scaled such that their sum is equal to L. If this field
%          is not provided or is entered as an empty matrix, then the links
%          will not be scaled.
%
%   physics: Structure containing information about the system's physics.
%       Fields are:
% 
%       drag_ratio: Ratio of lateral to longitudinal drag
%
%       drag_coefficient: Drag per unit length for longitudinal direction
%
%   jointangles: Angles of the joints between the links, defining the
%       chain's current shape.
%
%
% Output:
%
%   Mp: A matrix (m x m matrix, where m is the number of elements in
%       the joint angle input). This matrix encodes the resistance of the
%       viscous medium to motion of the joints. It acts as a linear map
%       from joint angular velocity to torques acting on the joints, and as
%       a quadratic map from joint angular velocity to power being
%       dissipated through the system. 
%
%       If Mp is taken as a metric tensor, pathlengths calculated from it
%       have lengths in units of "sqrt(power)*time". This pathlength is a
%       pacing-indepenent cost for the motion, and represents a "minimum
%       cost" for following the path, which is achieved for a trajectory
%       that follows the path at constant speed as measured by the metric,
%       which itself corresponds to constant power. Given two paths, the path
%       with a shorter length according to this metric can always be
%       traversed more quickly at a given power, or with a lower power in a
%       given time.


    %%%%%%%
    % 	To calculate this metric, first get the Local connection and full
    % 	mass matrix for the system

    [A, ~,~, ~,~,M_full] = Inertial_local_connection(geometry,physics,jointangles);
    
    % Now fold down the mass matrix onto the shape space
    
    Mp = [-A.' eye(size(A,2))] * M_full * [-A; eye(size(A,2))];

%     %%%%%%%%
%     % Now calculate the metric contribution to each link
%     
%     % Pre-allocate storage for the metric contributions. Each contribution
%     % is m x m (square matrix with as many rows/columns as there are shape
%     % variables), and there is one contribution per link.
%     link_metrics = repmat({zeros(numel(jointangles))},size(J_full));
% 
%     % Get the Metric contribution from each link. This is the drag matrix
%     % acting on an individual link, pulled back onto the space of joint
%     % angles by the total Jacobian (including locomotion effects) from
%     % joint angular velocity to the motion of the link.
%     for idx = 1:numel(link_metrics)
%         
%         link_metrics{idx} = LowRE_link_metric(J_full{idx},...            % Jacobian from body velocity of base link and shape velocity to body velocity of this link
%                                               A,...                      % Local connection (i.e. Jacobian from shape velocity to body velocity of base link)
%                                               h.lengths(idx),...         % Length of this link
%                                               physics.drag_ratio,...     % Ratio of lateral to longitudinal drag
%                                               physics.drag_coefficient); % Bulk drag coefficient
%     end
%     
%     %%%%%%%
%     % Finally, sum the metric contributions of the links to get the metric
%     % for the whole system.
%     Mp = sum(cat(3,link_metrics{:}),3);
% 
%     
%     
% end
% 
% function link_metric = LowRE_link_metric(J_full,A,L,drag_ratio,c)
% % Calculate the contribution to the power metric from a link. This is the
% % drag resistance on the link, pre- and post-multiplied by the Jacobian
% % from shape velocity to this link's body velocity
% 
% 	%%%%%%%
%     % Local drag on a link.
%     %   Longitudinal drag is proportional to longitudinal velocity. 
%     %   Lateral drag is proportional to lateral velocity, with a different drag
%     %       coefficient, specified by its ratio to the longitudinal drag. 
%     %   Rotational drag is proportional to rotational velocity, with a
%     %       coefficient based on integrating lateral drag along the length
%     %       of the link
%     %   c is the absolute drag coefficient
%     
%     % (Note that drag terms need to be positive here to make this a
%     % positive-definite matrix that we can use as a metric)
%      drag_matrix =     [L       0                0;
%                         0    drag_ratio*L       0;
%                         0        0           drag_ratio/12*L^3]*c;
%                     
%                     
%     %%%%%%%%%%%
%     % Calculate the jacobian from shape variables to body velocity of this
%     % point on the system:
%     
%     % -A maps shape velocity to body velocity of the system, so augmenting
%     % -A with an identity matrix produces a map from shape velocity to body
%     % and shape velocity of the system
%     J_intermediate = [-A; eye(size(A,2))];
%     
%     % J_full maps the system's body and shape velocity to the body velocity
%     % of this link.
%     J_total = J_full * J_intermediate;
%     
%     % Pre- and post-multiplying J_total by the drag matrix pulls the drag
%     % matrix back from being a metric on the link's body velocities to
%     % being a metric on the system's shape velocities
%     link_metric = J_total' * drag_matrix * J_total;
% 	
% end