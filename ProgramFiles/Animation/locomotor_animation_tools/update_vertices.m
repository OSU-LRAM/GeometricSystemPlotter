function object = update_vertices(object)
%Take an object and update the location of its global-position vertices
%based on its current global position and its local geometry.

%First check whether object has moved since the vertices were last 
% calcualted. If not, return the object untouched.
% Otherwise, store the current position to be used as "last position" in
% future calls.
% if ~isempty(object.ref_pos.last.calc)
%         
%     if (object.ref_pos.now == object.ref_pos.last.calc)
% 
%         return;
% 
%     end
%     
% else
%     
%     object.ref_pos.last.calc = object.ref_pos.now;
%     
% end


%ensure reference position is a row vector
if size(object.ref_pos.now,1)>size(object.ref_pos.now,2)
    object.ref_pos.now = object.ref_pos.now';
end

%split xy-position and theta from reference position
xy = object.ref_pos.now(1:2);
theta = object.ref_pos.now(3);

%rotation matrix
SE2 = vec_to_mat_SE2(object.ref_pos.now);%[cos(theta) -sin(theta); sin(theta) cos(theta)];

%Depending on object type, run the appropriate transform
switch object.type
    
%     case 'circle'
%  
%         %basis for circle parameterization
%         z = linspace(0,2*pi,30)';
% 
%         %unit circle primative
%         unit_circle = [cos(z) sin(z)];
%         unit_circle(end,:)=unit_circle(1,:);
% 
%         %Scale and translate circle
%         object.geom.vertices = (ones(size(z)) * xy)...
%             + (object.local_geom.radius * unit_circle);
%         
%     case 'rectangle'
%         
%         %ensure sidelengths is a column vector
%         if size(object.local_geom.sidelengths,1) < ...
%                 size(object.local_geom.sidelengths,2)
%             object.local_geom.sidelengths = object.local_geom.sidelengths';
%         end
%         
%         %unit square
%         unit_square = .5*[  1   1
%                             -1  1
%                             -1  -1
%                             1   -1
%                             1   1];
% 
% 
%         %rectangle at origin, unrotated
%         rectangle = unit_square * diag(object.local_geom.sidelengths);
% 
%         %rectangle at origin, rotated
%         rectangle = [SO2 * rectangle']';
% 
%         %rectangle translated, rotated
%         object.geom.vertices = (ones(size(unit_square)) * diag(xy))...
%             + rectangle;
        
    case 'polygon'

%         %ensure border wraps all the way around
%         if object.local_geom.vertices(1,1) ~=...
%                 object.local_geom.vertices(end,1)...
%                     ||...
%                         object.local_geom.vertices(1,2) ~=...
%                             object.local_geom.vertices(end,2)
%             
%             object.local_geom.vertices = ...
%                 [object.local_geom.vertices;object.local_geom.vertices(1,:)];
%         end

        %scale polygon
        polygon = object.local_geom.vertices;% * object.local_geom.scaling;

        %translate polygon
        if iscell(object.local_geom.vertices)
            for idx = 1:numel(object.local_geom.vertices)
                object.geom.vertices{idx} = (SE2 * polygon{idx}')';
            end
        else
            
            object.geom.vertices = (SE2 * polygon')';
        end
        
%     case 'line'
%         
%         %rotate and translate line
%         object.geom.vertices = ones(size(object.local_geom.vertices)) * ...
%             diag(xy) + object.local_geom.vertices;
        
    otherwise
        
        disp([object.type,' is not defined'])
        
end

% %find the maximum radius of the object, to prune collision detection.
% object.local_geom.max_rad = max(distance(object.ref_pos.now(1:2),object.geom.vertices));

if isfield(object.local_geom,'sense')
        
    %basis for circle parameterization
    z = linspace(0,2*pi,15)';

    %unit circle primative
    unit_circle = [cos(z) sin(z)];

    %Scale and translate circle
    object.geom.sense = (ones(size(z)) * xy)...
        + (object.local_geom.sense * unit_circle);
end