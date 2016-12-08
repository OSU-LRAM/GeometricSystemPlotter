function [B,h] = fatbackbone_from_links(linklengths,r,L,width,orientation)

    % Specify orientation as midpoint-tangent unless specified otherwise
    if ~exist('orientation','var')
        orientation = 'midpoint-tangent';
    end

    [h] = backbone_from_curvature_bases(kappa_basis_input,r,L);
    
    h_points = h(linspace(L*-.5,L*.5,100))';
	
	
	B = fatbackbone(h,L*[-.5 .5],width);
    
    % Rotate and translate locus if specified
    
    %Augment the backbone locus with 1s to make them SE(2) point vectors
    B_aug = [B,ones(size(B,1),1)];
    
    switch orientation
        case 'midpoint-tangent'
            % do nothing, locus is already in this condition
        case 'com-mean'
    
            % Build a rotation matrix for center of mass and mean orientation
            x = mean(h_points(:,1));
            y = mean(h_points(:,2));
            theta = mean(h_points(:,3));
            
            R = [cos(theta) -sin(theta) x;
                sin(theta) cos(theta) y;
                0 0 1];
            
            % Multiply inverse of this matrix by point location matrix
            
            B_aug = (R\B_aug')';
            
            
    end

    % Restore B to its un-augmented form
    B = B_aug(:,1:2);

end