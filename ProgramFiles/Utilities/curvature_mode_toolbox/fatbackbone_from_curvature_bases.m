function [B,h] = fatbackbone_from_curvature_bases(kappa_basis_input,r,L,width,orientation)

    % Specify orientation as midpoint-tangent unless specified otherwise
    if ~exist('orientation','var')
        orientation = 'midpoint-tangent';
    elseif or(isequal(orientation,'com-mean'), isequal(orientation,'midpoint-tangent'))
        % do nothing; orientation is a valid option
    else
        % Get the user folder path
        load('sysplotter_config.mat','inputpath');
        % The hypothetical file path, supposing orientation is a system name
        calcfilePath = fullfile(inputpath, 'sysplotter_data', ['sysf_', orientation, '_calc.mat']);
        if exist(calcfilePath, 'file')
            orientation = 'from-sys';
        else
            orientation = 'midpoint-tangent'; % illustrate_backbone_shapespace provides a warning about this
        end
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
            
        case 'from-sys'
		
            % Import the datafile of the specified system
            load(calcfilePath,'s');
            
            r_cell = num2cell(r);
            % find the value of Bx for the given shape variables:
            x = interpn(s.grid.eval{:}, s.B_optimized.eval.Beta{1}, r_cell{:});
            % find the value of By for the given shape variables:
            y = interpn(s.grid.eval{:}, s.B_optimized.eval.Beta{2}, r_cell{:});
            % find the value of Btheta for the given shape variables:
            theta = interpn(s.grid.eval{:}, s.B_optimized.eval.Beta{3}, r_cell{:});

            %%% Same as in com-mean:%%%
            R = [cos(theta) -sin(theta) x;
                sin(theta) cos(theta) y;
                0 0 1];
            % Multiply inverse of this matrix by point location matrix
            B_aug = (R\B_aug')';
    end

    % Restore B to its un-augmented form
    B = B_aug(:,1:2);

end
