function [h, J] = backbone(geometry,shapeparams)

% Generate backbone geometry and Jacobian from its local definition
switch geometry.type
    
    case {'curvature basis','curvature bases'}
        
        [h, J] = backbone_from_curvature_bases(geometry.function,shapeparams,geometry.length);
        
    case 'general curvature'
        
        [h, J] = backbone_from_general_curvature(geometry.function,shapeparams,geometry.length);
        
    otherwise
        warning('backbone type not supported')
        
end

end