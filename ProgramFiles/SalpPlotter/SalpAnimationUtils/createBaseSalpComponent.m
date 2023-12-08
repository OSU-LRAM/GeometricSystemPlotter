function comp = createBaseSalpComponent(salp, comp_color, structSize)
% Create structure to store salp components

    if nargin < 3
        structSize = salp.geometry.nLinks;
    end
    
    comp = struct();
    comp.color = comp_color;
    comp.points = cell(structSize, 1);
end