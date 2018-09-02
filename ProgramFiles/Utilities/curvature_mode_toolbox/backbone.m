function [h, J] = backbone(geometry,shapeparams)

% If no baseframe is specified, use a centered chain
if ~isfield(geometry,'baseframe') || isempty(geometry.length)
    baseframe = 'centered';
else
    baseframe = geometry.baseframe;
end

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