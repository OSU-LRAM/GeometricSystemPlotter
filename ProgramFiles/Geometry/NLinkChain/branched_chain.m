function [h, J, J_full,frame_zero,J_zero] = branched_chain(geometry,shapeparams_global)
% Build a backbone for a branching chain of links. The geometry input
% should be a cell array of structures, each of which has the same elements
% as required for an N_link_chain, and also (except for the first chain)
% has an "attachment" field that specifies the chain it should be attached
% to and a location at which to attach it, using the same format as
% baseframe.
%
% Inputs:
%
%   geometry: Structure containing information about geometry of chain.
%       Fields are:
%
%       linklengths: Vectors of link lengths defining a kinematic chain.
%           For a subchain (any chain other than the first), setting the
%           first link length to zero will mean that it directly hinges off
%           of the chain to which it is attached, rather than having a
%           fixed link attached to the chain
%
%       baseframe (optional): Specification of which link on the body should be
%           treated as the base frame. Default behavior is to put the base
%           frame at the center of the primary chain, either at the midpoint of the
%           middle link, or at the joint between the two middle links and at
%           the mean of their orientation. Options for this field are:
% 
%               'centered' :    Default behavior
%               'tail' :        Lowest-numbered link is the base frame
%               'tail-tip' :    Start of lowest-numbered link is the base frame
%               'head' :        Highest-numbered link is the base frame
%               'head-tip' :    End of the highest-numbered link is the base frame
%               numeric :       Specify a link number to use as a base frame
%               sysf_           Pull minimum-perturbation coordinates from a
%                                       sysf_ file. Argument should be the name of
%                                       a system in the current UserFiles folder
%               'start' :       Modifier on a link specification (e.g., 
%                   {2,'start'} to put the frame at the proximal end of the
%                   specified link
%               'end' :         Modifier on a link specification (e.g., 
%                   {head,'end'} to put the frame at the distal end of the
%                   specified link
%               transform:      A 3x3 SE(2) matrix giving the position
%                   of the base frame. This can be a standalone entry or a
%                   modifier on any other frame.
%
%       modes (optional): Option to map input "shapeparams" across
%           multiple links in a chain (which can have more joints than the
%           provided number of "shapeparams". Should be a matrix in which
%           each column is the set of joint angles on the full chain associated
%           with a unit value of the corresponding "jointangle" input.
%
%
%       length (optional): Total length of the chain. If specified, the elements of
%           will be scaled such that their sum is equal to L. If this field
%           is not provided or is entered as an empty matrix, then the links
%           will not be scaled.
%
%   shapeparams: A vector of the angles between the links. Must be one
%       element shorter than the linklengths vector%
%   
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
%           mat_to_vec_SE2
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
%   frame_zero : The transformation from the first link to the selected
%           baseframe
%
%   J_zero : The Jacobian from shape velocity to body velocity of the
%       selected baseframe, with the first link held fixed.
%       

% Build a lookup table from joint numbers to joints on the subchains
joint_lookup_local_to_global = cell(size(geometry.subchains));
joint_count = 0;

joint_lookup_global_to_local = [];

for idx = 1:numel(joint_lookup_local_to_global)
    
    % Handle any local modes applied to the system
    % If no modes are specified, use an identity mapping for the modes
    if ~isfield(geometry.subchains{idx},'modes') || isempty(geometry.subchains{idx}.modes)
        geometry.subchains{idx}.modes = eye(numel(geometry.subchains{idx}.linklengths)-1);
    end

    joint_lookup_local_to_global{idx} = joint_count + (1:size(geometry.subchains{idx}.modes,2))';
    
    joint_count = joint_count + size(geometry.subchains{idx}.modes,2);
    
    joint_lookup_global_to_local = [joint_lookup_global_to_local;                                                           % Anything already in joint_lookup
                    [size(joint_lookup_global_to_local,1)+(1:size(geometry.subchains{idx}.modes,2))' ...     % A sequence of numbers starting where existing numbers leave off
                , idx*ones(size(geometry.subchains{idx}.modes,2),1) ...                        % The subchain number
                , (1:size(geometry.subchains{idx}.modes,2))'] ];                             % The joint number within the subchain
                 
    
end

% Handle any global modes applied to the system
% If no modes are specified, use an identity mapping for the modes
if ~isfield(geometry,'modes') || isempty(geometry.modes)
    modes = eye(numel(shapeparams_global));
else
    modes = geometry.modes;
end

% Expand jointangles from specified shape variables to actual joint angles
% by multiplying in the modal function
shapeparams_global = shapeparams_global(:);
shapeparams = modes*shapeparams_global;


% Get the kinematics for each individual chain
[h_set, J_set, J_full_set,frame_zero_set,J_zero_set,chain_description_set] ...
    = cellfun(@(geom,shape) N_link_chain(geom,shapeparams(shape)),geometry.subchains,joint_lookup_local_to_global,'UniformOutput',false);

%%%%%%%%%
% Internally pad the Jacobians with zeros so that they line up with their
% actual joint numbers

J_set_padded = cell(size(J_set));
J_full_set_padded = cell(size(J_full_set));
J_temp_set_padded = cell(size(chain_description_set));

% Iterate over the individual chains
for idx = 1:numel(J_set)
    
    
    % Extract the local-to-global joint lookup for joint on this subchain 
    jllg = joint_lookup_local_to_global{idx};
    
    
    % Create matrices of the right size to hold the J and J full matrices
    J_set_padded{idx} = repmat({zeros(3,numel(shapeparams))},size(J_set{idx}));
    J_full_set_padded{idx} = repmat({zeros(3,3+numel(shapeparams))},size(J_set{idx}));
    J_temp_set_padded{idx} = repmat({zeros(3,numel(shapeparams))},size(chain_description_set{idx}.J_temp));
        
    % Iterate over the Jacobans in the chain
    for idx2 = 1:numel(J_set_padded{idx})
        
        % Move the columns of J_set to their global joint number
        J_set_padded{idx}{idx2}(:,jllg) = J_set{idx}{idx2};
    
        % Move the shape columns of the J_full_set
        J_full_set_padded{idx}{idx2}(:,1:3) = J_full_set{idx}{idx2}(:,1:3);
        J_full_set_padded{idx}{idx2}(:,3+jllg) = J_full_set{idx}{idx2}(:,4:end);
         
        % Move the columns of J_full_set to their global joint number
        J_temp_set_padded{idx}{idx2}(:,jllg) = chain_description_set{idx}.J_temp{idx2};

    end
   
    % Overwrite the original J values with their padded versions
    J_set{idx} = J_set_padded{idx};
    J_full_set{idx} = J_full_set_padded{idx};
    chain_description_set{idx}.J_temp = J_temp_set_padded{idx};
    
end

% Loop over the subchains
for idx = 2:numel(h_set)

    % Extract the attachment parameters for ths subchain, and put it into
    % the chain_description
	attach = geometry.subchains{idx}.attach;
    
    chain_description_parent = chain_description_set{attach.parent};
    
    chain_description_parent.baseframe = attach.location;
    
    % Put the padded J_temp into the chain_description
    chain_description_set{idx}.J_temp = J_temp_set_padded{idx};
    
    % Use the modified chain descriptions to get the location and Jacobian of the
    % attachment point relative to the first chain's baseframe
    [frame_zero,J_zero,link_zero] = N_link_conversion_factors(chain_description_parent);
    
    % Move the chain to frame_zero, J_zero
    [h_m,J_set{idx},J_full_set{idx},chain_description_set{idx}] ...
        = N_link_conversion_move_chain(chain_description_set{idx},frame_zero,J_zero);
    
    h_set{idx}.pos = mat_to_vec_SE2(h_m);
    
    % Update the cumulative sum of joint angles
    if ~iscell(attach.location)
        attach.location = {attach.location};
    end
    for idx2 = 1:numel(attach.location)
        
        if isscalar(attach.location{idx2})
            orientation_offset = chain_description_parent.jointangles_c(attach.location{idx2});
        elseif ischar(attach.location{idx2}) && strncmp(attach.location{idx2},'tail',4)
            orientation_offset = chain_description_parent.jointangles_c(1);
        elseif ischar(attach.location{idx2}) && strncmp(attach.location{idx2},'head',4)
            orientation_offset = chain_description_parent.jointangles_c(end);
        elseif isnumeric(attach.location{idx2})
            t_v = mat_to_vec_SE2(attach.location{idx2});
            orientation_offset = t_v(3);
        elseif ischar(attach.location{idx2}) && ( strncmp(attach.location{idx2},'start',4) || strncmp(attach.location{idx2},'end',4))
            orientation_offset = 0;
        else
            orientation_offset = 0;
            warning('Chain attachment does not call out a link, errors may occur if the whole chain is specified to be in averaged coordinates')
        end
        chain_description_set{idx}.jointangles_c ...
            = chain_description_set{idx}.jointangles_c ...
             + orientation_offset;
    end

    
%     % Extract the attachment parameters for ths subchain
% 	attach = geometry{idx}.attach;
%     
% 	% Get the position of the link to which this subchain is attached
% 	h_attach = vec_to_mat_SE2(h_set{attach.parent}.pos(attach.link,:));
%     
%     % Multiply the link transformations by the attachment transformation
%     for idx2 = 1:size(h_set{idx}.pos,2)
%         h_set{idx}.pos(idx2,:) = mat_to_vec_SE2( h_attach * vec_to_mat_SE2(h_set{idx}.pos(idx2,:)));
%     end
%     
%     % Transform the Jacobians
%     for idx2 = 1:numel(J_set_padded{idx})
%         
%         % Use left lifted action of the attachment point to modify the
%         % in-frame Jacobian
%         J_set_padded{idx}{idx2} = TeLg(h_attach) * J_set_padded{idx}{idx2};
%         
%         % Convert the full Jacobian
%         Adjointinverse_body_transform = Adjinv(h_attach); % Adjoint-inverse transformation by position of attachment point
%         
%         J_full_set_padded{idx}{idx2} = [Adjointinverse_body_transform * J_full_set_padded{idx}{idx2}(:,1:3), J_full_set_padded{idx}{idx2}(:,4:end)];
%         
%     end
    
end


% % Combine the h_set values into a single item
% h.pos = [];
% h.lengths = [];
% for idx = 1:numel(h_set)
%     h.pos = [h.pos;h_set{idx}.pos];
%     h.lengths = [h.lengths;h_set{idx}.lengths];
% end
% 
% J = {};
% for idx = 1:numel(J_set)
%     
%     J = [J,J_set{idx}]; %#ok<AGROW>
%     
% end
% 
% 
% J_full = {};
% for idx = 1:numel(J_set)
%     
%     J_full = [J_full,J_full_set{idx}]; %#ok<AGROW>
%     
% end


%%%%%%%%%%%
% Apply any global baseframe and mode transformations



% Baseframe modifications
chain_description.chain_m = [];
chain_description.jointchain_m = [];
chain_description.links_m = [];
chain_description.joints_m = [];
chain_description.links_v = [];
chain_description.joints_v = [];
chain_description.jointangles = [];
chain_description.jointangles_c = [];
chain_description.linklengths = [];
chain_description.shapeparams = [];
chain_description.modes = [];
chain_description.J_temp = [];
chain_description.baseframe = {};

for idx = 1:numel(chain_description_set)
    
    chain_description.chain_m = cat(3,chain_description.chain_m, chain_description_set{idx}.chain_m);
    chain_description.jointchain_m = cat(3,chain_description.jointchain_m, chain_description_set{idx}.jointchain_m);
    chain_description.links_m = cat(3,chain_description.links_m, chain_description_set{idx}.links_m);
    chain_description.joints_m = cat(3,chain_description.joints_m, chain_description_set{idx}.joints_m);
    chain_description.links_v = [chain_description.links_v; chain_description_set{idx}.links_v];
    chain_description.joints_v = [chain_description.joints_v; chain_description_set{idx}.joints_v];
    chain_description.jointangles = [chain_description.jointangles; chain_description_set{idx}.jointangles];
    chain_description.jointangles_c = [chain_description.jointangles_c; chain_description_set{idx}.jointangles_c];
    chain_description.linklengths = [chain_description.linklengths; chain_description_set{idx}.linklengths];
    chain_description.J_temp = [chain_description.J_temp, chain_description_set{idx}.J_temp];
    
end
    
chain_description.shapeparams = shapeparams;
chain_description.modes = modes;
chain_description.baseframe = geometry.baseframe;

% Make the final conversion by the overall baseframe specification
[frame_zero,J_zero] = N_link_conversion_factors(chain_description);
    
[h_m,J,J_full,chain_description] = N_link_conversion(chain_description,frame_zero,J_zero); 


% Apply modal restriction to shape components of Jacobians
for idx = 1:numel(J)
    J{idx} = J{idx} * modes;
    J_full{idx}(:,4:end) = J_full{idx}(:,4:end)*modes;
end


% For output, convert h into row form. Save this into a structure, with
% link lengths included
h.pos = mat_to_vec_SE2(h_m);

h.lengths = [];
for idx = 1:numel(h_set)
    h.lengths = [h.lengths;h_set{idx}.lengths];
end


end