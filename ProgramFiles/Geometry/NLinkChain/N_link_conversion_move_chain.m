function [h_m,J,J_full,C] = N_link_conversion_move_chain(C,frame_zero,J_zero)
%%%%%%%
% This is a helper-function for N_link_chain. 
%
% This function moves the existing baseframe of the chain described in C to
% a frame at frame_zero relative to the baseframe. (This is different from
% N_link_conversion, which leaves the chain fixed and moves the baseframe
%
% N_link_chain takes in an argument that specifies which link on the chain
% (or another frame, such as a center-of-mass frame or one loaded from a
% system file).
%
% All of the N_link_chain calculations are made with the lowest-numbered
% (leftmost, or 'tail') link as the base link.
%
% Transforming the positions of the link relative to the base frame is
% straightforward: We simply multiply each transformation by the inverse of
% the transformation from the old base frame to the new base frame
%
% To transform the link Jacobians into the new coordinates, we need the
% Jacobian from joint velocities to body velocity of new frame, with
% first link fixed. This calculation depends on how the kinematic map
% from the original base frame to the new base frame is calculated.
%
% The code here uses the conversion transformation and Jacobian calculated
% in N_link_conversion_factors.


    % Extract elements from structure array (some functions like size do not
    % work properly if applied to struct contents without a breakout)
    chain_m = C.chain_m;
    jointchain_m = C.jointchain_m;
    J_temp = C.J_temp;

    % Preallocate a matrix for holding the transformed link and joint location matrices
    h_m = zeros(size(chain_m));
    hj_m = zeros(size(jointchain_m));
    
    % If we're working with symbolic variables, then we need to explicitly make
    % the array symbolic, because matlab tries to cast items being inserted
    % into an array into the array class, rather than converting the array to
    % accomodate the class of the items being inserted 
    if isa(chain_m,'sym')
        h_m = sym(h_m);
    end

    if isa(jointchain_m,'sym')
        hj_m = sym(hj_m);
    end
    
    
    % For each link and joint in the chain, transform its matrix by the inverse of the
    % central transformation
    for idx = 1:size(h_m,3)
        h_m(:,:,idx) = frame_zero * chain_m(:,:,idx);
        
        if isa(chain_m,'sym')
            h_m(:,:,idx) = simplify(h_m(:,:,idx),'steps',10);
        end
        
    end
    
    for idx = 1:size(hj_m,3)
        hj_m(:,:,idx) = frame_zero * jointchain_m(:,:,idx);
        
        if isa(chain_m,'sym')
            hj_m(:,:,idx) = simplify(hj_m(:,:,idx),'steps',10);
        end
        
    end

    
    % Save h_m to the chain description structure as the new chain_m
    C.chain_m = h_m;
    C.jointchain_m = hj_m;


    %%%%%%%%
    % Use the Jacobian for the transformation from the old base frame to
    % the new base frame to compute the link Jacobians with respect to that
    % frame
    
    % First operation is to get the Jacobian that maps shape velocity to
    % body velocity of the links, taking the frame that frame_zero is
    % defined with respect to as fixed
    %
    % This is achieved by multiplying the Jacobian from the original base
    % frame to the new base frame by the Adjoint-inverse of the
    % transformation from the new base frame to the link, and then
    % subtracting this value from the original 
    J_new = J_temp;
    for idx = 1:numel(J_new)
        J_new{idx} = ...
            (J_temp{idx} + ...                % Jacobian from joint velocities to body velocity of link, with first link of the chain fixed
                Adjinv(chain_m(:,:,idx)) * ...       % Adjoint-inverse transformation by position of this link relative to frame_zero
                    J_zero);                 % Jacobian from joint velocities to body velocity of frame_zero, with its reference frame fixed
    end
    
    % Save this Jacobian out to the chain description structure
    C.J_temp = J_new;


    % Second operation is to transform the link Jacobians such that they
    % return the velocities of the links in the coordinate directions of
    % the new base frame (the reference frame for frame_zero), instead of in each link's local coordinates
    %
    % This is achieved by multiplying each link Jacobian by the left lifted
    % action of the transformation from the new base frame to the link
    J = J_new;
    for idx = 1:numel(J_new)
        J{idx} = TeLg(h_m(:,:,idx)) * J_new{idx}; % Left lifted action rotates into new base frame coordinates
    end


    % Third operation is to calculate the full Jacobian from body velocity
    % of the base frame and shape velocity to body velocity of each link.
    %
    % This is achieved by augmenting each of the link's body-velocity
    % Jacobians with a new block that contains the Adjoint-inverse of the
    % transformation from the new base frame to the link
    J_full = J_new;
    for idx = 1:numel(J_new)

        Adjointinverse_body_transform = Adjinv(h_m(:,:,idx)); % Adjoint-inverse transformation by position of this link relative to frame zero

        J_full{idx} = [Adjointinverse_body_transform J_full{idx}]; % Place the Adjoint-inverse matrix as the first three columns of J_full
    end

    

end