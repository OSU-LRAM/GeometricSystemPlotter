function [h,J,J_full] = backbone_conversion(h_orig,J_orig,frame_zero,J_zero)
%%%%%%%
% This is a helper-function for backbone. 
%
% backbone takes in an argument that specifies which link on the chain
% (or another frame, such as a center-of-mass frame or one loaded from a
% system file).
%
% All of the backbone calculations are made with the midpoint
% (s=0) frame as the base frame.
%
% Transforming the positions of the backbone pieces relative to the base
% frame is straightforward: We simply multiply the transformation function
% for the backbone pieces by the inverse of the transformation from the old
% base frame to the new base frame
%
% To transform the link Jacobians into the new coordinates, we need the
% Jacobian from shape velocities to body velocity of new frame, with
% the center frame fixed. This calculation depends on how the kinematic map
% from the original base frame to the new base frame is calculated.
%
% The code in this file handles the special cases required for finding the
% transformation and Jacobian from the old base frame to the new base frame
% for different specifications of the new base frame.
%
% The code here uses the conversion transformation and Jacobian calculated
% in backbone_conversion_factors.

% Only calculate Jacobians if requested in the output
if nargout > 1
    calc_J = 1;
else
    calc_J = 0;
end


    % New h function premultiplies the original h function by the
    % transformation to the new basis
    h = @(s)backbone_conversion_helper(s,h_orig,frame_zero);
    

    if calc_J
        
        % First operation is to get the Jacobian that maps shape velocity to
        % velocity of the links, taking our new base frame as fixed. 
        %
        % This is achieved by multiplying the Jacobian from the original base
        % frame to the new base frame by the Adjoint-inverse of the
        % transformation from the new base frame to the link, and then
        % subtracting this value from the original. We then multiply these
        % quantities by TeLg(frame_zero) to convert from body velocity of
        % the links to their world velocity

        J = @(s) TeLg(frame_zero) * (J_orig(s) - Adjinv(h(s))*J_zero);

        % Second operation is to calculate the full Jacobian from body velocity
        % of the base frame and shape velocity to body velocity of each link.
        %
        % This is achieved by augmenting each of the link's body-velocity
        % Jacobians with a new block that contains the Adjoint-inverse of the
        % transformation from the new base frame to the link

        J_full = @(s) [Adjinv(h(s)) (J_orig(s) - Adjinv(h(s))*J_zero)];
        
    end
    

end


function h_out = backbone_conversion_helper(s,h_orig,frame_zero)

    % Evaluate original backbone function at each s value provided and
    % convert to SE(2) matrix form
    h_m = vec_to_mat_SE2(h_orig(s)');
    
    % Loop over the sheets in the stack, multiplying each by the inverse of
    % frame_zero
    for idx = 1:size(h_m,3)
        h_m(:,:,idx) = frame_zero \ h_m(:,:,idx);
    end
    
    % Return the new frame locations as columns of an array
    h_out = mat_to_vec_SE2(h_m)';


end