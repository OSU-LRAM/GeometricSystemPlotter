function [h, J, J_full] = backbone_from_links(linklengths,jointangles,L,baseframe)
% Build a backbone for a chain of links, specified as a vector of link
% lengths and the joint angles between them.
%
% Inputs:
%
%   linklengths: A vector of link lengths defining a kinematic chain
%
%   jointangles: A vector of the angles between the links. Must be one
%       element shorter than the linklengths vector
%
%   L (optional): Total length of the chain. If specified, the elements of
%       will be scaled such that their sum is equal to L. If this argument
%       is not provided or is entered as an empty matrix, then the links
%       will not be scaled.
%   
%   baseframe (optional): Specification of which link on the body should be
%       treated as the base frame. Default behavior is to put the base
%       frame at the center of the chain, either at the midpoint of the
%       middle link, or at the joint between the two middle links and at
%       the mean of their orientation. Options for this argument are:
%       
%           'centered' :    Default behavior
%           'tail' :        Lowest-numbered link is the base frame
%           'head' :        Highest-numbered link is the base frame
%           'head-tip':     End of the highest-numbered link is the base frame
%           numeric :       Specify a link number to use as a base frame
%
%
% Outputs:
%
%   h : Locations of the chain links relative to the selected base frame,
%           and lengths of the links
%
%       h.pos: The link locations are stored in an Nx3 array, with
%           each row a link and the columns corresponding to x,y,theta
%           values for that link. These can be converted to stacks of SE(2)
%           matrices via vec_to_mat_SE2, and back to vectors by
%           mat_to_vec_SE(2)
%   
%   	h.lengths : Vector of linklengths. This is the original link length
%           specification passed into the function, scaled by L, and stored
%           as a column of values.
%
%   J : Jacobians from joint angle velocities to velocity of each link
%           relative to  the base frame (the standard derivative of h).
%           These Jacobians are stored in a cell array, with one matrix per
%           cell.
%
%   J_full : Full Jacobians from body velocity of the base frame and joint
%           angle velocities to body velocity of each link relative to a
%           non-moving frame.
%

%%%%%%%%%%%%
% Input parsing

% If no length is specified, do not scale the links
if ~exist('L','var') || isempty(L)
    L = sum(linklengths);
end

% If no baseframe is specified, use a centered chain
if ~exist('baseframe','var')
    baseframe = 'centered';
end

% Decide if we need to calculate the Jacobian
calc_J = nargout > 2;
calc_J_full = nargout >3;

% Force linklength and jointangle vectors to be columns, and normalize link
% lengths for a total length of 1.
linklengths = linklengths(:)/sum(linklengths)*L;
jointangles = jointangles(:);

%%%%%%%%%%%%



%%%%%%%%%%%%%%
% Get some basic information about the chain whose kinematics we're
% calculating

% Number of links and joints
N_links = numel(linklengths);
M_joints = numel(jointangles);
    if M_joints ~= N_links-1
        error(['There should be one more link than joint. There are ' num2str(N_links) ' links and ' num2str(M_joints) ' joint angles specified']);
    end
% Decide if there is an even or odd number of joints
N_odd = mod(N_links,2);

%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%
% Core kinematics code

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

%%%%
% Starting with the first link, work along the chain to find the location
% of each midpoint

% Create a 3-d array in which the ith 2-dimensional sheet is the 3x3 SE(2) matrix
% representing the corresponding link's configuration
chain_m = repmat(eye(3),1,1,N_links);

% The first link is taken as the identity for now, so we iterate over each
% of the following links
for idx = 2:N_links

    % The transformation from one midpoint to the next is the product
    % of the half-link transformation on the proximal link, the joint
    % transformation, and the half-link transformation on the distal
    % link. 

    chain_m(:,:,idx) = chain_m(:,:,idx-1) * ... % Each link's position is based off the position of the previous link
        links_m(:,:,idx-1)*...                  % Transform along the proximal link
        joints_m(:,:,idx-1)*...                 % Rotate by the intermediate joint angle
        links_m(:,:,idx);                       % Transform along the distal link 

end


%%%%%%%%
% If we will be calculating the Jacobian or have an even number of links,
% also calculate the location of each joint
if (calc_J) || (~N_odd)
    
    % Build a 3-d array to hold the SE(2) matrices for each joint
    % (specifically for the frame at the end of the link proximal to the
    % joint, treating that as the stator for the joint).
    jointchain_m = repmat(eye(3),1,1,numel(jointangles));
    
    % Place the first joint at the end of the first link
    jointchain_m(:,:,1) = links_m(:,:,1);
    
    % Iterate along the chain to find the configuration of the later joints
    for idx = 2:numel(jointangles)

        % The transformation from one joint to the next is the product
        % of the proximal joint's transformation with the link between the two joints 

        jointchain_m(:,:,idx) = ...
            jointchain_m(:,:,idx-1) * ... % Each joint's position is based off the position of the previous link
            joints_m(:,:,idx-1) * ...      % Rotate by the angle of the previous joint
            links_m(:,:,idx)^2;           % Move along the link twice (because our transforms are half-links)

    end

end

%%%%%%%%%%%%%%
%Code for finding the Jacobian

% Only calculate the Jacobian if requested
if calc_J

    % Initially, we calculate the Jacobians from shape velocity to body
    % velocity of each link, with the first link held fixed. We later
    % transform this into other Jacobians for the system that are more
    % interesting
    
    % Calculate a Jacobian for each link
    J_temp = repmat({zeros(3,M_joints)},1,N_links);
    for idx = 1:N_links

        % For each Jacobian, calculate one column for each joint
        for idx2 = 1:M_joints

            % Links are only sensitive to the motion of joints proximal to
            % them in the chain, so check if the link index is higher than
            % the joint index
            sensitivity = idx > idx2; 
            
            if sensitivity

                % Calculate the displacement of the link relative to the
                % joint as the inverse of the current joint's location,
                % multiplied by the current link's location
                relative_transform = jointchain_m(:,:,idx2) \ chain_m(:,:,idx);

                % Find the Adjoint-inverse transformation corresponding to
                % this relative position
                Adjointinverse_transform = Adjinv(relative_transform);
                

                % Multiply these the Adjoint-inverse transformation by the
                % joint axis to get the Jacobian from the current joint to
                % the current link
                J_temp{idx}(:,idx2) = Adjointinverse_transform * [0;0;1];

            else

                % If this link is not sensitive to this joint, leave this
                % column of its Jacobian as zero
                
            end

        end
        

    end
    
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For final output, convert link locations to be referenced off of a
% desired base frame

    %%%%%%%
    % To transform the link Jacobians into the new coordinates, we need the
    % Jacobian from joint velocities to body velocity of new frame, with
    % first link fixed. This calculation depends on how the kinematic map
    % from the original base frame to the new base frame is calculated. 
    %
    % In particular, if we're putting the frame on the chain, we can use
    % the same iteration-down-the-chain algorithm we used above. If we're
    % putting the frame at a floating average of the link positions, we can
    % combine the Jacobians for the individual links into a Jacobian for
    % the averaged frame


% Depending on which baseframe option is specified, use different methods
% to calculate transform and Jacobian to use in the conversion
switch baseframe
    
    % Places the reference frame at the lowest-index link on the chain
    case 'tail'
        
        % Identify the first link
        link_zero = 1;
        
        % Extract the transform from the end link to itself 
        frame_zero = eye(3);       
        
        % Jacobian for base link is zero
        J_zero = zeros(size(J_temp{1}));

    % Places the reference frame at the middle of the chain, splitting the
    % difference between the two middle links if there is an even number of
    % links
    case {'centered','center'}


        % If there is an odd number of links, the middle link is the center of the
        % chain
        if N_odd

            % Identify the middle link
            link_zero = ceil(N_links/2);

            % Extract the transform from the end link to the middle link
            frame_zero = chain_m(:,:,link_zero);
            
            %%%%%%
            % Jacobian for conversion is Jacobian of link zero
            J_zero = J_temp{link_zero};

        % If there is an even number of links, the midpoint is at the middle joint,
        % rotated by half of that joint's angle
        else

            % Identify the middle joint
            joint_zero = ceil(M_joints/2);

            % Extract the position of the middle joint, and multiply it by half of
            % the transform associated with its joint angle
            frame_zero = jointchain_m(:,:,joint_zero) * ...
                                vec_to_mat_SE2(joints_v(joint_zero,:)/2);
                            
            %%%%%%%%
            % Jacobian for conversion is Jacobian of link before it, with a
            % half-step to get to the end of the link, and a half rotation
            % to average orientation between the two links
            halfstep = Adjinv(links_v(joint_zero,:));
            halfrotation = Adjinv(joints_v(joint_zero,:)/2);
            J_zero = halfrotation * halfstep * J_temp{joint_zero};
            
            J_zero(:,joint_zero) = [0; 0; .5]; % Account for centered frame being half-sensitive to middle joint
           

        end
        
        
%         %%%%%%%%%%%%%%%%
%         % Calculate columns of Jacobian to new base frame 
%         for idx2 = 1:M_joints
% 
%             % New frame is  only sensitive to the motion of joints proximal to it in the chain
%             if N_odd
%                sensitivity = link_zero > idx2;  
%             else
% 
%                 % If the frame is proximal to the joint, it is not sensitive to it
%                 if idx2 > joint_zero
%                     sensitivity = 0;
%                 % If the frame is at the center joint, it has half-sensitivity to the joint
%                 elseif idx2 == joint_zero 
%                     sensitivity = 0.5;
%                 % If the frame is distal to the joint, it is sensitive to it
%                 else
%                     sensitivity = 1;
%                 end
% 
%             end
% 
% 
%             if sensitivity
% 
%                 % Calculate the displacement of the new frame relative to the joint
%                 relative_transform = jointchain_m(:,:,idx2)\frame_zero;
% 
%                 % Calculate the adjoint-inverse transformation corresponding to
%                 % the joint-to-link transformation
%                 Adjointinverse_transform = Adjinv(relative_transform);
% 
% 
%                 % Multiply these transformations by the joint axis to get the
%                 % Jacobian column (multiply by sensitivity to catch case of
%                 % half-sensitivity
%                 J_zero(:,idx2) = sensitivity * Adjointinverse_transform * [0;0;1];
% 
%             else
% 
%                 % If this link is not sensitive to this joint, give it a
%                 % column of zeros in its Jacobian
%                 J_zero(:,idx2) = zeros(3,1);
% 
%             end
% 
%         end
        
        
        
    % Places the reference frame on the highest-numbered link in the chain
    case 'head'
        
        % Identify the end link
        link_zero = N_links;

        % Extract the transform from the end link to the end link
        frame_zero = chain_m(:,:,link_zero);
        
        %%%%%%%%
        % Jacobian to new frame is Jacobian of last link
        J_zero = J_temp{end};

    % Places the reference frame at the *end of* the highest-numbered link
    % in the chain
    case 'head-tip'
        
        % Identify the end link
        link_zero = N_links;

        % Extract the transform from the end link to the end link, and
        % multiply it by the half-link transformation on that link
        frame_zero = chain_m(:,:,link_zero)*links_m(:,:,link_zero);
        
        %%%%%%%%
        % Jacobian to new frame is Jacobian of last link, but with an
        % adjoint-inverse transform by the half-link to get to the end
        halfstep = Adjinv(links_v(end,:));
        J_zero = halfstep * J_temp{end};

    % Places the reference frame at center of mass and average orientation
    % of the links, using link-lengths as weighting terms
    case 'com-mean'
        
        % Convert link positions to row form
        chain = mat_to_vec_SE2(chain_m);
        
        % Take a weighted average of the link positions
        CoM = sum(diag(linklengths)*chain)/L;
        
        % Place the new frame at this location
        frame_zero = vec_to_mat_SE2(CoM);

        %%%%%%%%%%%
        % The Jacobian of the weighted average of frames is the
        % weighted average of their Jacobian (by the commutativity of
        % sumation and derivation operations).

        % Multiply each link's Jacobian by its link length
        J_weighted = J_temp;
        for idx = 1:numel(J_weighted)
            J_weighted{idx} = TeLg(chain(idx,:)) * J_temp{idx} * linklengths(idx);
        end

        % Sum the weighted Jacobians
        J_zero = sum(cat(3,J_weighted{:}),3);

        % Divide by the total length to g
        J_zero = J_zero/L;  
        
        % Bring into local coordinates
        J_zero = TgLginv(frame_zero)*J_zero;
        
        
    otherwise
        
        % Check for numeric baseframe specification
        if isnumeric(baseframe)
                        
            % Make sure that base frame specification is actually in
            % the range of valid link numbers
            if baseframe <= N_links &&  baseframe > 0

                % Set the base frame as numerically specified in the input
                link_zero = baseframe;

                % Extract the transform from the end link to the end link, and
                % multiply it by the half-link transformation on that link
                frame_zero = chain_m(:,:,link_zero);

                %%%%%%%%
                % Jacobian to new frame is Jacobian of nth link
                J_zero = J_temp{baseframe};

            else

                error (['Baseframe specification ' num2str(baseframe) ' is not in the range of link numbers'])

            end

                
            
        else 
       
            error (['Baseframe specification ' baseframe ' is not a valid option'])
            
        end
        
end
        
        
        
        

% Preallocate a matrix for holding the transformed link matrices
h_m = zeros(size(chain_m));

% For each link in the chain, transform its matrix by the inverse of the
% central transformation
for idx = 1:size(h_m,3)
    h_m(:,:,idx) = frame_zero \ chain_m(:,:,idx);
end

% For output, convert h into row form. Save this into a structure, with
% link lengths included
h.pos = mat_to_vec_SE2(h_m);
h.lengths = linklengths;

%

%%%%%%%%%%%%



if calc_J
    
    %%%%%%%
    % To transform the link Jacobians into the new coordinates, we need the
    % Jacobian from joint velocities to body velocity of new frame, with
    % first link fixed. This calculation depends on how the kinematic map
    % from the original base frame to the new base frame is calculated. 
    %
    % In particular, if we're putting the frame on the chain, we can use
    % the same iteration-down-the-chain algorithm we used above. If we're
    % putting the frame at a floating average of the link positions, we can
    % combine the Jacobians for the individual links into a Jacobian for
    % the averaged frame

    % Decide which method to use for calculating the Jacobian from the old
    % frame to the new frame
%         switch new_frame_jacobian
% 
%             % Build the Jacobian for the new frame by working along the chain
%             % from link 1 to the specified location
%             case 'on chain'
% 
%                 % Calculate column of Jacobian with respect to each joint proximal
%                 % to the new base frame
%                 for idx2 = 1:M_joints
% 
%                     % New frame is  only sensitive to the motion of joints proximal to it in the chain
%                     if N_odd
%                        sensitivity = link_zero > idx2;  
%                     else
% 
%                         % If the frame is proximal to the joint, it is not sensitive to it
%                         if idx2 > joint_zero
%                             sensitivity = 0;
%                         % If the frame is at the center joint, it has half-sensitivity to the joint
%                         elseif idx2 == joint_zero 
%                             sensitivity = 0.5;
%                         % If the frame is distal to the joint, it is sensitive to it
%                         else
%                             sensitivity = 1;
%                         end
% 
%                     end
% 
% 
%                     if sensitivity
% 
%                         % Calculate the displacement of the new frame relative to the joint
%                         relative_transform = jointchain_m(:,:,idx2)\frame_zero;
% 
%                         % Calculate the adjoint-inverse transformation corresponding to
%                         % the joint-to-link transformation
%                         Adjointinverse_transform = Adjinv(relative_transform);
% 
% 
%                         % Multiply these transformations by the joint axis to get the
%                         % Jacobian column (multiply by sensitivity to catch case of
%                         % half-sensitivity
%                         J_zero(:,idx2) = sensitivity * Adjointinverse_transform * [0;0;1];
% 
%                     else
% 
%                         % If this link is not sensitive to this joint, give it a
%                         % column of zeros in its Jacobian
%                         J_zero(:,idx2) = zeros(3,1);
% 
%                     end
% 
%                 end
% 
%             case 'averaged'
% 
%                 % The Jacobian of the weighted average of frames is the
%                 % weighted average of their Jacobian (by the commutativity of
%                 % sumation and derivation operations).
% 
%                 % Multiply each link's Jacobian by its link length
%                 J_weighted = J_temp;
%                 for idx = 1:numel(J_weighted)
%                     J_weighted{idx} = J_temp{idx} * linklengths(idx);
%                 end
% 
%                 % Sum the weighted Jacobians
%                 J_zero = sum(cat(3,J_weighted{:}),3);
% 
%                 % Divide by the total length to g
%                 J_zero = J_zero/L;
% 
%         end

    %%%%%%%%
    % Now that we have the Jacobian for the new base frame, we can compute
    % the link Jacobians with respect to that frame
    
    % First operation is to get the Jacobian that maps shape velocity to
    % body velocity of the links, taking our new base frame as fixed. 
    %
    % This is achieved by multiplying the Jacobian from the original base
    % frame to the new base frame by the Adjoint-inverse of the
    % transformation from the new base frame to the link, and then
    % subtracting this value from the original 
    J_new = J_temp;
    for idx = 1:numel(J_new)
        J_new{idx} = ...
            (J_new{idx} - ...                % Jacobian from joint velocities to body velocity of link, with first link fixed
                Adjinv(h.pos(idx,:)) * ...       % Adjoint-inverse transformation by position of this link in new frame
                    J_zero);                 % Jacobian from joint velocities to body velocity of new frame, with first link fixed
    end


    % Second operation is to transform the link Jacobians such that they
    % return the velocities of the links in the coordinate directions of
    % the new base frame, instead of in each link's local coordinates
    %
    % This is achieved by multiplying each link Jacobian by the left lifted
    % action of the transformation from the new base frame to the link
    J = J_new;
    for idx = 1:numel(J_new)
        J{idx} = TeLg(h.pos(idx,:)) * J_new{idx}; % Left lifted action rotates into new base frame coordinates
    end


    % Third operation is to calculate the full Jacobian from body velocity
    % of the base frame and shape velocity to body velocity of each link.
    %
    % This is achieved by augmenting each of the link's body-velocity
    % Jacobians with a new block that contains the Adjoint-inverse of the
    % transformation from the new base frame to the link
    J_full = J_new;
    if nargout > 2 % Only make this calculation if J_full is requested

        for idx = 1:numel(J_new)

            Adjointinverse_body_transform = Adjinv(h.pos(idx,:)); % Adjoint-inverse transformation by position of this link relative to frame zero

            J_full{idx} = [Adjointinverse_body_transform J_full{idx}]; % Place the Adjoint-inverse matrix as the first three columns of J_full
        end

    end

end


end