function [B,h_fun,J] = fatbackbone_from_stretch_curvature(curvdef,cparams,L,width,orientation)

    % Specify orientation as midpoint-tangent unless specified otherwise
    if ~exist('orientation','var')
        orientation = 'midpoint-tangent';
    end
    
    if nargout == 3
        [h,J] = backbone_from_stretchable_curvature(curvdef,cparams,L);
    else  %if nargout == 2
        h = backbone_from_stretchable_curvature(curvdef,cparams,L);
    end
    % 	    h_points = h(linspace(L*-.5,L*.5,100))';  
        
% 	B = fatbackbone(h,L*[-.5 .5],width);

    % Rotate and translate locus if specified
    
    %Augment the backbone locus with 1s to make them SE(2) point vectors
%   B_aug = [B,ones(size(B,1),1)];
    
    switch orientation
        case 'midpoint-tangent'
            % do nothing, locus is already in this condition
            
            B = fatbackbone(h,L*[-.5 .5],width);
            h_fun = @(s) h(s);
        case 'com-mean'

            h_points = h(linspace(L*-.5,L*.5,100))';
            
            % Build a rotation matrix for center of mass and mean orientation
            x = mean(h_points(:,1));
            y = mean(h_points(:,2));
            theta = mean(h_points(:,3));
            
            R = [cos(theta) -sin(theta) x;
                sin(theta) cos(theta) y;
                0 0 1];
            
            B = fatbackbone(h,L*[-.5 .5],width);
            % augment B and h
            B_aug = [B,ones(size(B,1),1)];
%             h_points_aug = [h_points(:,1:2),ones(size(h_points,1),1)];
            
            % Multiply inverse of this matrix by point location matrix

            h_fun = {@(s) h(s);
                     @(h) [cos(mean(h(:,3))) -sin(mean(h(:,3))) mean(h(:,1));...
                           sin(mean(h(:,3)))  cos(mean(h(:,3))) mean(h(:,2));...
                           0 0 1]};
            B_aug = (R\B_aug')';
            B = B_aug(:,1:2);
            
            % Jacquelin's code implements like below when called 
            % pts is the 1xn array containing points along the backbone to be
            %   plotted
            % h_fun outputs two function handles, the first gets the
            %   location on the backbone and the second gets the rotation
            %   matrix for that location
            % h_arrows is my final output array for the x and y points to plot
            % h_arrows_rot calls the second function to get the right
            %   rotation matrix given where on the backbone we are
            % h_aug gets the current backbone positions of each point in
            %   the pts array
            
%             h_arrows_rot = h_fun{2}(h_fun{1}(linspace(L*-0.5,L*0.5,100))');
%             h_aug = h_fun{1}(pts);
%             h_arrows = (h_arrows_rot\[h_aug(1:2,:);ones(1,size(h_aug,2))])';
            
    end

    % Restore B to its un-augmented form
%     B = B_aug(:,1:2);

end