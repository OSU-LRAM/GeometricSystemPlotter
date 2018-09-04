function [B,h] = fatbackbone_from_general_curvature(curvdef,cparams,L,width,orientation)

% This is a legacy function, rewritten to be a wrapper on the new
% general-purpose fat backbone

geometry.type = 'general curvature';
geometry.function = curvdef;
geometry.length = L;
geometry.baseframe = orientation;

display.aspect_ratio = width/L;

shapeparams = cparams;

[B,h] = fat_backbone(geometry,shapeparams,display);

%     % Specify orientation as midpoint-tangent unless specified otherwise
%     if ~exist('orientation','var')
%         orientation = 'midpoint-tangent';
%     end
%     
%     [h] = backbone_from_general_curvature(curvdef,cparams,L);
%     h_points = h(linspace(-.5,.5,100))';
% 	
% 	B = fatbackbone(h,[-.5 .5],width);
% 
%     % Rotate and translate locus if specified
%     
%     %Augment the backbone locus with 1s to make them SE(2) point vectors
%     B_aug = [B,ones(size(B,1),1)];
%     
%     switch orientation
%         case 'midpoint-tangent'
%             % do nothing, locus is already in this condition
%         case 'com-mean'
%     
%             % Build a rotation matrix for center of mass and mean orientation
%             x = mean(h_points(:,1));
%             y = mean(h_points(:,2));
%             theta = mean(h_points(:,3));
%             
%             R = [cos(theta) -sin(theta) x;
%                 sin(theta) cos(theta) y;
%                 0 0 1];
%             
%             % Multiply inverse of this matrix by point location matrix
%             
%             B_aug = (R\B_aug')';
%             
%         otherwise
%             
%             % Make the string input robust against the user possibly
%             % including sysf_ in the function
%             if startsWith(orientation,'sysf_')
%                 orientation = orientation(6:end);
%             end
%             
%             load('sysplotter_config.mat','inputpath')
%             calcfilePath = fullfile(inputpath, 'sysplotter_data',...
%                 ['sysf_', orientation, '_calc.mat']);
%             
%             if exist(calcfilePath, 'file') % Confirm that the system file exists
%                 % Import the datafile of the specified system
%                 load(calcfilePath,'s')
% 
%                 r_cell = num2cell(cparams);
%                 % find the value of Bx for the given shape variables:
%                 x = interpn(s.grid.eval{:}, s.B_optimized.eval.Beta{1},r_cell{:});
%                 % find the value of By for the given shape variables:
%                 y = interpn(s.grid.eval{:}, s.B_optimized.eval.Beta{2}, r_cell{:});
%                 % find the value of Btheta for the given shape variables:
%                 theta = interpn(s.grid.eval{:}, s.B_optimized.eval.Beta{3}, r_cell{:});
% 
%                 %%% Same as in com-mean:%%%
%                 R = [cos(theta) -sin(theta) x;
%                     sin(theta) cos(theta) y;
%                     0 0 1];
%                 % Multiply inverse of this matrix by point location matrix
%                 B_aug = (R\B_aug')';
%                 
%             else % If the file doesn't exist, show a warning
%                 warning(strcat('system data file sysf_',orientation,...
%                     '_calc.mat can not be found. ',...
%                     'Defaulting to midpoint-tangent orientation.'))
%                 % then do nothing, to match midpoint-tangent case above.
%             end
%     end
% 
%     % Restore B to its un-augmented form
%     B = B_aug(:,1:2);

end
