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
        calcfilePath = strcat(inputpath,'\sysplotter_data\sysf_',orientation,'_calc.mat');
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
            
            % Calculate tolerance for find() (this will vary depending on
                                            % system evaluation density)
            tol=abs(s.grid.eval{1}(1,1)-s.grid.eval{1}(2,1))*0.5;
            
            % Determine the location of the specified values of a1 and a2
            [k,j]=find(abs(s.grid.eval{1}-r(1))<tol); % find indicies where a1 is as specified
            row=k(1); % because a1 is represented in the x axis, k will be a vector of identical values. save just one
            [k,j]=find(abs(s.grid.eval{2}(1,j)-r(2))<tol); % find the indicies where a1 AND a2 are as specified
            col=j(1);
            x = s.B_optimized.eval.Beta{1}(row,col); % find the value of Bx for the given a1 and a2
            y = s.B_optimized.eval.Beta{2}(row,col); % find the value of By for the given a1 and a2
            theta = s.B_optimized.eval.Beta{3}(row,col); % find the value of Btheta for the given a1 and a2
            
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