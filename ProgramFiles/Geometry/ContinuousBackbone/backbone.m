function [h, J, J_full,frame_zero,J_zero] = backbone(geometry,shapeparams)
% Build a continuous backbone from a curvature function
%
% Inputs:
%
%   geometry: Structure containing information about geometry of chain.
%       Fields used by this function are:
%
%       type: The kind of backbone function provided here. Used by backbone
%           to calculate the geometry of the backbone
%
%       function: A handle to a curvd_ function describing the system's
%           curvature function, or a cell vector of function handles to the
%           curvature modes
%
%       baseframe (optional): Specification of which link on the body should be
%           treated as the base frame. Default behavior is to put the base
%           frame at the center of the chain, either at the midpoint of the
%           middle link, or at the joint between the two middle links and at
%           the mean of their orientation. Options for this field are:
% 
%               'centered' :    Default behavior
%               'tail' :        Base frame at s = -0.5
%               'head' :        Base frame at s =  0.5
%               numeric :       Specify a position on the range [-0.5 0.5]
%                                   to use as the base frame
%               sysf_           Pull minimum-perturbation coordinates from a
%                                       sysf_ file. Argument should be the name of
%                                       a system in the current UserFiles folder
%
%       length (optional): Total length of the backbone. If specified, the elements of
%           will be scaled such that their sum is equal to L. If this field
%           is not provided or is entered as an empty matrix, then the
%           backbone will not be scaled.
%
%   shapeparams: A vector of the  input parameters taken by the curvature
%       function
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



% If no baseframe is specified, use a centered backbone
if ~isfield(geometry,'baseframe') || isempty(geometry.baseframe)
    baseframe = 'centered';
else
    baseframe = geometry.baseframe;
end

% If no length is specified, use unit length
if ~isfield(geometry,'length') || isempty(geometry.length)
    geometry.length = 1;
end


% Only calculate Jacobians if requested in the output
if nargout > 1
    calc_J = 1;
else
    calc_J = 0;
end

% Generate backbone geometry and Jacobian from its local definition
switch geometry.type

    case {'curvature basis','curvature bases'}

        if calc_J
            [h, J] = backbone_from_curvature_bases(geometry.function,shapeparams,geometry.length);
        else
            h = backbone_from_curvature_bases(geometry.function,shapeparams,geometry.length);
        end

    case 'general curvature'

        if calc_J
            [h, J] = backbone_from_general_curvature(geometry.function,shapeparams,geometry.length);
        else
            h = backbone_from_general_curvature(geometry.function,shapeparams,geometry.length);
        end
    otherwise
        warning('backbone type not supported')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Code above this block is kinematics of a backbone. Code below this block is
% for changing the base frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        
%%%%%%
% Calculate the transformation from the original base frame to the new base
% frame

if calc_J
    [frame_zero,J_zero] = backbone_conversion_factors(h,J,shapeparams,baseframe,geometry.length);   
else
    frame_zero = backbone_conversion_factors(h,[],shapeparams,baseframe,geometry.length);
end

%%%%%%
% Use frame_zero and J_zero to convert the link transformations and
% Jacobian so that they are refefenced off of the new base frame
if calc_J
    [h,J,J_full] =backbone_conversion(h,J,frame_zero,J_zero);
else
    h = backbone_conversion(h,[],frame_zero,[]); 
end


end