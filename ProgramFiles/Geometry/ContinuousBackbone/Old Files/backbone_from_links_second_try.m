function [h, J, J_full] = backbone_from_links_second_try(linklengths,jointangles,L)
% Build a backbone for a chain of links.

% If no length is specified, do not scale the links
if ~exist('L','var')
    L = sum(linklengths);
end

% Number of links and joints
N_links = numel(linklengths);
M_joints = numel(jointangles);
    if M_joints ~= N_links-1
        error(['There should be one more link than joint. There are ' num2str(N_links) ' links and ' num2str(M_joints) ' joint angles specified']);
    end
% Decide if there is an even or odd number of joints
N_odd = mod(N_links,2);
    
    
% Force linklength and jointangle vectors to be columns, and normalize link
% lengths for a total length of 1.
linklengths = linklengths(:)/sum(linklengths)*L;
jointangles = jointangles(:);

%%%%%
% Convert the links and joints to SE(2) matrices

% Divide the link lengths by two, since we're going midpoint to
% midpoint
halflengths = linklengths/2;
    
% Convert the link half-length values to [l 0 0] rows, encoding group
% elements displaced along the x direction from the identity.
links_v = [halflengths, zeros(N_links,2)];


% Convert the joint angle values to group elements displaced along the
% theta direction from the identity
joints_v = [zeros(M_joints,2), jointangles];

% Find the matrix representations of the link transformations
links_m = vec_to_mat_SE2(links_v);

% Find the matrix representations of the joint transformations
joints_m = vec_to_mat_SE2(joints_v);


% Now work along the chain to find the location of each midpoint
chain_m = repmat(eye(3),1,1,N_links);
for idx = 2:N_links

    % The transformation from one midpoint to the next is the product
    % of the half-link transformation on the proximal link, the joint
    % transformation, and the half-link transformation on the distal
    % link. 

    chain_m(:,:,idx) = chain_m(:,:,idx-1) * ...                   %Each link's position is based off the position of the previous link
        links_m(:,:,idx-1)*joints_m(:,:,idx-1)*links_m(:,:,idx);  % half-joint-half product

end


%%%%%%%%
% If we will be calculating the Jacobian or have an even number of links,
% also calculate the location of each joint
if (~N_odd) || (nargout>1)
    
    jointchain_m = repmat(eye(3),1,1,numel(jointangles));
    jointchain_m(:,:,1) = links_m(:,:,1);
    for idx = 2:numel(jointangles)

        % The transformation from one joint to the next is the product
        % of the proximal joint's transformation with the link between the two joints 

        jointchain_m(:,:,idx) = jointchain_m(:,:,idx-1) * ... %Each link's position is based off the position of the previous link
            joints_m(:,:,idx-1)*links_m(:,:,idx)^2;  % joint-half-half product

    end

end

%%%%%
% For final output, convert into middle-link coordinates


% If there is an odd number of links, the middle link is the center of the
% chain
if N_odd
    
    % Identify the middle link
    link_zero = ceil(N_links/2);
    
    % Extract the transform from the end link to the middle link
    h_transform = chain_m(:,:,link_zero);
    
% If there is an even number of links, the midpoint is at the middle joint,
% rotated by half of that joint's angle
else
    
    % Identify the middle joint
    joint_zero = ceil(M_joints/2);
    
    % Extract the position of the middle joint, and multiply it by half of
    % the transform associated with its joint angle
    h_transform = jointchain_m(:,:,joint_zero)*vec_to_mat_SE2(joints_v(joint_zero,:)/2);

end

% Preallocate a matrix for holding the transformed link matrices
h_m = zeros(size(chain_m));

% For each link in the chain, transform its matrix by the inverse of the
% central transformation
for idx = 1:size(h_m,3)
    h_m(:,:,idx) = h_transform \ chain_m(:,:,idx);
end

% For output, convert h into row form.
h = mat_to_vec_SE2(h_m);

%%%%%%%%%%%%%%
%Code for finding the Jacobian

% Only calculate the Jacobian if requested
if nargout > 1

    % Transform joint locations to centered frame
    for idx = 1:size(jointchain_m,3)
        jointchain_mc(:,:,idx) = h_transform \ jointchain_m(:,:,idx);
    end

    % Calculate Jacobian for each link
    for idx = 1:N_links

        % Calculate column of Jacobian with respect to each joint
        for idx2 = 1:M_joints

            sensitivity = J_influence(idx2,idx,M_joints,N_links);
            if sensitivity

                % Calculate the displacement of the link relative to the joint
                relative_transform = jointchain_mc(:,:,idx2)\h_m(:,:,idx);

                % Calculate the adjoint-inverse transformation corresponding to
                % the joint-to-link transformation
                Adjointinverse_transform = Adjinv(mat_to_vec_SE2(relative_transform));
                
                % Calculate the lifted left transformation to rotate
                % coordinates from local-link to base-link
                Leftlifted_transform = TeLg(mat_to_vec_SE2(h_m(:,:,idx)));

                % Multiply these transformations by the joint axis to get the
                % Jacobian column
                J{idx}(:,idx2) = Leftlifted_transform ...
                    * Adjointinverse_transform * sensitivity * [0;0;1];

                % If a third output is requested, it is the Jacobian from
                % system body and shape velocity to body velocity of the
                % joint
                if nargout > 2
                    

                    % Multiply the Adjoint transformation by the joint axis to get the
                    % Jacobian column
                    J_full{idx}(:,idx2) = Adjointinverse_transform * sensitivity * [0;0;1];

                end
            else

                % If this link is not sensitive to this joint, give it a
                % column of zeros in its Jacobian
                J{idx}(:,idx2) = zeros(3,1);
                
                % Likewise for the full Jacobian
                if nargout > 2
                    J{idx}(:,idx2) = zeros(3,1);
                end

            end

        end
        
        % If calculating the full Jacobian, put the block for motion of the
        % whole system here
        if nargout > 2
            
            Adjointinverse_from_center = Adjinv(mat_to_vec_SE2(h_m(:,:,idx)));
            Leftlifted_inverse_transform = TgLginv(mat_to_vec_SE2(h_m(:,:,idx)));
            
            J_full{idx} = [Adjointinverse_from_center (Leftlifted_inverse_transform * J{idx})];
            
        end

    end
end
    

end


function sensitivity = J_influence(m,n,M_joints,N_links)
% Return the sensitivity of link n to joint m when the middle of the chain
% is taken as fixed.

    % Behavior if there is an odd number of links
    if mod(N_links,2)
        
        % Identify middle link
        link_zero = ceil(N_links/2);
        joint_separator = (M_joints+1)/2;
        
        % If the link is further back than the joint, it rotates negatively
        % when the joint rotates positively. Since we're fixing the middle
        % link as our body frame, this condition only applies for joints to
        % the left of the middle link.
        if (n <= m) && (m < joint_separator)
            
            sensitivity = -1;
            
        % If the link is further along than the joint, it rotates
        % positively when the joint rotates positively. Since we're fixing
        % the middle link as our body frame, this condition only applies
        % for links to the right of the middle link.
        elseif n > m && (m > joint_separator)
            
            sensitivity = 1;
            
        % If the link is closer to the middle link than the joint, it is
        % not affected by motion of the link
        else
            
            sensitivity = 0;
            
        end
        
    % Behavior if there is an even number of links    
    else
        
        joint_zero = ceil(M_joints/2);
        joint_separator = (M_joints+1)/2;
        
        % If the link is further back than the joint, it rotates negatively
        % when the joint rotates positively. Since we're fixing the middle
        % link as our body frame, this condition only applies for links to
        % the left of the middle joint.
        if (n <= m) && (m <= joint_separator)
            
            sensitivity = -1;
            
        % If the link is further along than the joint, it rotates
        % positively when the joint rotates positively. Since we're fixing
        % the middle link as our body frame, this condition only applies
        % for links to the right of the middle link.
        elseif n > m && (m >= joint_separator)
            
            sensitivity = 1;
            
        % If the link is closer to the middle link than the joint, it is
        % not affected by motion of the link
        else
            
            sensitivity = 0;
            
        end
        
        % If the joint is the middle joint, multiply the sensitvity by 1/2,
        % since its influence is split across the two halves of the system.
        if m == joint_zero
            sensitivity = 0.5*sensitivity;
        end
        
    end

end

