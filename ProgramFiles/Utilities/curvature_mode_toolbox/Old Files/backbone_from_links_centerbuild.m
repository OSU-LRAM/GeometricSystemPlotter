function [h, J, J_full] = backbone_from_links_centerbuild(linklengths,jointangles,L)
% Build a backbone for a chain of links.


% Force linklength and jointangle vectors to be columns, and normalize link
% lengths for a total length of 1.
linklengths = linklengths(:);%/sum(linklengths);
jointangles = jointangles(:);

%%%%%
% Split the chain into two halves, at either the middle of the middle link,
% or at the middle joint (respectively for odd and even numbers of links

% Determine number of links and decide if it is odd or even
N_links = numel(linklengths);
N_odd = mod(N_links,2) ~= 0;

% Number of links in each side chain
n_links = floor(N_links/2);

% Links in chains to the left and right of the center
links{1} = linklengths(n_links:-1:1);
links{2} = linklengths(end-n_links+1:end);

% If the number of links is odd, the middle link is link zero. Otherwise,
% add a zero-length middle link;
if N_odd
    link_zero = linklengths(n_links+1);
else
    link_zero = 0;
end

% Add link zero to the base of each chain
links{1} = [link_zero;links{1}];
links{2} = [link_zero;links{2}];


%%%%%
% Now divide up the joint angles

% Calculate number of joints and compare it against the number of links
M_joints = numel(jointangles);
    if M_joints ~= N_links-1
        error(['There should be one more link than joint. There are ' num2str(N_links) ' links and ' num2str(M_joints) ' joint angles specified']);
    end
M_odd = mod(M_joints,2) ~= 0 ;


% Number of joints in each side chain
m_joints = floor(M_joints/2);

% Links in chains to the left and right of the center
joints{1} = jointangles(m_joints:-1:1);
joints{2} = jointangles(end-m_joints+1:end);

% If the number of joints is odd, split the middle joint across the two
% side chains
if M_odd 
    middle_angle = jointangles(m_joints+1);
    joints{1} = [middle_angle/2; joints{1}];
    joints{2} = [middle_angle/2; joints{2}];
end

%%%%%
% Now flip the signs on the left chain (because we're working right to left)
links{1} = -links{1};
joints{1} = -joints{1};

%%%%%
% Now convert the links and joints to SE(2) matrices

%Iterate over the left and right subchains
for idx = 1:numel(links)
    
    % Convert the link length values to group elements displaced along the
    % x direction from the identity.
    %
    % Divide the link lengths by two, since we're going midpoint to
    % midpoint
    links_v{idx} = [links{idx}/2, zeros(numel(links{idx}),2)];
    
    
    % Convert the joint angle values to group elements displaced along the
    % theta direction from the identity
    joints_v{idx} = [zeros(numel(joints{idx}),2), joints{idx}];
    
    
    % Find the matrix representations of the link transformations
    links_m = vec_to_mat_SE2(links_v{idx});
    
    % Find the matrix representations of the joint transformations
    joints_m = vec_to_mat_SE2(joints_v{idx});
    
    
    % Now work along the chain to find the location of each midpoint
    chain_m{idx} = repmat(eye(3),1,1,numel(links{idx}));
    for idx2 = 2:numel(links{idx})
    
        % The transformation from one midpoint to the next is the product
        % of the half-link transformation on the proximal link, the joint
        % transformation, and the half-link transformation on the distal
        % link. 
        
        chain_m{idx}(:,:,idx2) = chain_m{idx}(:,:,idx2-1) * ...                   %Each link's position is based off the position of the previous link
            links_m(:,:,idx2-1)*joints_m(:,:,idx2-1)*links_m(:,:,idx2);  % half-joint-half product
                
    end
    
    % Make the chain values into rows of group element values instead of
    % representations
    chain{idx} = mat_to_vec_SE2(chain_m{idx}); %#ok<AGROW>
    
    % Remove the first element of the chain structure -- it is either not
    % actually a link (even numbers of links base off of a joint) or
    % appears on both chains, and will be added back in when gluing them
    % together.
    chain{idx} = chain{idx}(2:end,:); %#ok<AGROW>
    
    
    %%%%%%%%
    % If we will be calculating the Jacobian, also calculate the location
    % of each joint
    jointchain_m{idx} = repmat(eye(3),1,1,numel(joints{idx}));
    jointchain_m{idx}(:,:,1) = links_m(:,:,1);
    for idx2 = 2:numel(joints{idx})
    
        % The transformation from one joint to the next is the product
        % of the proximal joint's transformation with the link between the two joints 
        
        jointchain_m{idx}(:,:,idx2) = jointchain_m{idx}(:,:,idx2-1) * ... %Each link's position is based off the position of the previous link
            joints_m(:,:,idx2-1)*links_m(:,:,idx2).^2;  % joint-half-half product
                
    end
    
    
end


%%%%%%%%%%%%%
% Now glue the chain pieces together

% If there is an odd number of links, add a link centered at zero.
if N_odd
    chain_zero = [0 0 0];
else
    chain_zero = [];
end

% Concatenate the chains along with a possible zero-link. Flip the ordering
% of the left chain while doing this, to return it to left-to-right
% ordering.
h = [chain{1}(end:-1:1,:);
     chain_zero;
     chain{2}];


%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%

%Code for finding the Jacobian
%Iterate over the left and right subchains
for idx = 1:numel(links)
    
    % Calculate Jacobian for each link in subchain, excepting the base link
    for idx2 = 2:numel(links{idx})
        
        % Calculate column of Jacobian for each joint in subchain
        for idx3 = 1:numel(joints{idx})
            
            % Calculate the displacement of the link relative to the joint
            relative_transform = jointchain_m{idx}(:,:,idx3)\chain_m{idx}(:,:,idx2);
            
            % Calculate the adjoint-inverse transformation corresponding to
            % this transformation
            Adjointinverse_transform = Adjinv(mat_to_vec_SE2(relative_transform));
            
            % Multiply this transformation by the joint axis to get the
            % Jacobian column
            J_half{idx}{idx2}(:,idx3) = Adjointinverse_transform * joints_v{idx}(idx3,:)';
            
        end
        
    end
    
    % Trim out the empty element for the base link Jacobian
    J_half{idx} = J_half{idx}(2:end);

end

% Augment Jacobians by prepending or postpending columns of zeros
% corresponding to joints on other chains

for idx = 1:numel(links)
    J_filler{idx} = zeros(size(J_half{idx}{1}));
    
    % Remove one column from J_filler if there are an odd number of joints
    if M_odd
        J_filler{idx} = J_filler{idx}(:,2:end);
    end
end

% iterate over the left and right subchains (for a multi-chain system,
% this would expand to iterating over all subchains)
for idx = 1:numel(links)
    
    % In each subchain corresponding to higher-indexed joints, prepend columns
    % of zeros corresponding to joints in this subchain
    for idx2 = idx+1:numel(links)
        
        % Prepend zeros to each Jacobian in the subchain being operated on
        for idx3 = 1:numel(J_half{idx2})
        
            J_half{idx2}{idx3} = [J_filler{idx}, J_half{idx2}{idx3}];
            
        end
        
    end

    % In each subchain corresponding to lower-indexed joints, postpend columns
    % of zeros corresponding to joints in this subchain
    for idx2 = idx-1:-1:1
        
        % postpend zeros to each Jacobian in the subchain being operated on
        for idx3 = 1:numel(J_half{idx2})
        
            J_half{idx2}{idx3} = [J_half{idx2}{idx3}, J_filler{idx}];
            
        end
        
    end

end

%%%%
% Glue the pieces of the Jacobian together

% Concatenate the chains along with a possible zero-link. Flip the ordering
% of the left chain while doing this, to return it to left-to-right
% ordering.
%
% If there is an odd number of links, add a link centered at zero.
if N_odd
    J_zero = zeros(size(J_half{1}{1}));
    
    J = [J_half{1}(end:-1:1,:),...
     J_zero,...
     J_half{2}];
    
% If there is an even number of links overlay the midjoints and divide them
% by 2 (because they are really one joint, and each chain only sees half
% the motion
else
    J = [J_half{1}(end:-1:2),...
     0.5*(J_half{1}{1}+J_half{2}{1}),...
     J_half{2}(2:end)];
end



 
 
end