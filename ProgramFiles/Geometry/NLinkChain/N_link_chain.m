function [h, J, J_full,frame_zero,J_zero,chain_description] = N_link_chain(geometry,shapeparams)
% Build a backbone for a chain of links, specified as a vector of link
% lengths and the joint angles between them.
%
% Inputs:
%
%   geometry: Structure containing information about geometry of chain.
%       Fields are:
%
%       linklengths: A vector of link lengths defining a kinematic chain
%
%       baseframe (optional): Specification of which link on the body should be
%           treated as the base frame. Default behavior is to put the base
%           frame at the center of the chain, either at the midpoint of the
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
%       modes (optional): Option to map input "jointangles" across
%           multiple links in a chain (which can have more joints than the
%           provided number of "jointangles". Should be a matrix in which
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
%       element shorter than the linklengths vector
%
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
%   frame_zero : The transformation from the first link to the selected
%           baseframe
%
%   J_zero : The Jacobian from shape velocity to body velocity of the
%       selected baseframe, with the first link held fixed.
%
%   chain_description : The full set of chain features passed to
%       N_link_conversion_factors (but modified so that they are in the
%       specified baseframe, not the first-link frame.
%       

%%%%%%%%%%%%
% Input parsing

% If no length is specified, do not scale the links
if ~isfield(geometry,'length') || isempty(geometry.length)
    L = sum(geometry.linklengths);
else
    L = geometry.length;
end

% If no baseframe is specified, use a centered chain
if ~isfield(geometry,'baseframe') || isempty(geometry.baseframe)
    baseframe = 'centered';
else
    baseframe = geometry.baseframe;
end

% If no modes are specified, use an identity mapping for the modes
if ~isfield(geometry,'modes') || isempty(geometry.modes)
    modes = eye(numel(shapeparams));
else
    modes = geometry.modes;
end

% Force linklength and shapeparam vectors to be columns, and normalize link
% lengths for a total length of 1.
linklengths = geometry.linklengths(:)/sum(geometry.linklengths)*L;
shapeparams = shapeparams(:);


% Expand jointangles from specified shape variables to actual joint angles
% by multiplying in the modal function
jointangles = modes*shapeparams;

%%%%%%%%%%%%

% Prevent Matlab from playing tricks with imaginary numbers on symbolic
% inputs and from complaining about assumptions on constants
if isa(jointangles,'sym')
        
    assume(symvar(jointangles),'real');

    warning('off','symbolic:sym:sym:AssumptionsOnConstantsIgnored')
end

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

% Also take a cumulative sum of the joint angles -- this is useful if we
% want to take a mean orientation of the link orientations while making a
% distinction between orientations separated by 2pi radians
jointangles_c = [0; cumsum(jointangles)];

% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if or( isa(jointangles,'sym'), isa(linklengths,'sym') )
    chain_m = sym(chain_m);
end

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
    
    % Simplify trigonometric expressions if this is being calculated
    % symbolically
    if isa(chain_m,'sym')
        chain_m(:,:,idx) = simplify(chain_m(:,:,idx),'steps',10);
    end

end


%%%%
% Starting with the first joint, work along the chain to find the location
% of each joint

% Build a 3-d array to hold the SE(2) matrices for each joint
% (specifically for the frame at the end of the link proximal to the
% joint, treating that as the stator for the joint).
jointchain_m = repmat(eye(3),1,1,numel(jointangles));

% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if or( isa(jointangles,'sym'), isa(linklengths,'sym') )
    jointchain_m = sym(jointchain_m);
end

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

    % Simplify trigonometric expressions if this is being calculated
    % symbolically
    if isa(chain_m,'sym')
        jointchain_m(:,:,idx) = simplify(jointchain_m(:,:,idx),'steps',10);
    end
    
end



%%%%%%%%%%%%%%
%Code for finding the Jacobian


% Initially, we calculate the Jacobians from shape velocity to body
% velocity of each link, with the first link held fixed. We later
% transform this into other Jacobians for the system that are more
% interesting

% Calculate a Jacobian for each link

% Array of the size of one Jacobian
J_pattern = zeros(3,M_joints);

% If we're working with symbolic variables, then we need to explicitly make
% the array symbolic, because matlab tries to cast items being inserted
% into an array into the array class, rather than converting the array to
% accomodate the class of the items being inserted 
if or( isa(jointangles,'sym'),isa(linklengths,'sym') )
    J_pattern = sym(J_pattern);
end

J_temp = repmat({J_pattern},1,N_links);
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

            % Prevent symbolic chain expression from getting to cumbersome
            if isa(Adjointinverse_transform,'sym')
                Adjointinverse_transform = simplify(Adjointinverse_transform);
            end


            % Multiply these the Adjoint-inverse transformation by the
            % joint axis to get the Jacobian from the current joint to
            % the current link
            J_temp{idx}(:,idx2) = Adjointinverse_transform * [0;0;1];

        else

            % If this link is not sensitive to this joint, leave this
            % column of its Jacobian as zero

        end

    end
    
    % Convert joint-angle Jacobian into shape-mode coordinates
    J_temp{idx} = J_temp{idx} * modes;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code above this block is kinematics of a chain. Code below this block is
% for changing the base link

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
%%%%%%
% Calculate the transformation from the original base frame to the new base
% frame

chain_description = struct( ...
'chain_m',{chain_m},...
'jointchain_m',{jointchain_m},...
'links_m',{links_m},...
'joints_m',{joints_m},...
'links_v',{links_v},...            % Future refactoring could dispense with passing both links_v and links_m and joints_v and joints_m
'joints_v',{joints_v},...
'jointangles',{jointangles},...
'jointangles_c',{jointangles_c}, ...
'linklengths',{linklengths},...
'shapeparams',{shapeparams},...
'modes',{modes},...
'J_temp',{J_temp},...
'baseframe',{baseframe} ...
);

[frame_zero,J_zero] = N_link_conversion_factors(chain_description);        

%%%%%%
% Use frame_zero and J_zero to convert the link transformations and
% Jacobian so that they are refefenced off of the new base frame
[h_m,J,J_full,chain_description] = N_link_conversion(chain_description,frame_zero,J_zero); 




% %%%%%%%
% % Multiply the Jacobians by the modal matrices to produce Jacobians that
% % act from the modal coefficients rather than the full joint space
% J = cellfun(@(j) j*modes,J,'UniformOutput',false);
% full_mode_conversion = [eye(size(J{1},1)), zeros(size(J{1},1),size(modes,2));
%                         zeros(size(modes,1),size(J{1},1)),modes];
% J_full = cellfun(@(j) j*full_mode_conversion,J_full,'UniformOutput',false);

% For output, convert h into row form. Save this into a structure, with
% link lengths included
h.pos = mat_to_vec_SE2(h_m);
h.lengths = linklengths;


if nargout > 5
    % Get the partial derivative of the Jacobian
    dJdq = mobile_jacobian_derivative(J_full);
end


end