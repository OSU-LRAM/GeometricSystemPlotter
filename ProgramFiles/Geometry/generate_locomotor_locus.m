function B = generate_locomotor_locus(geometry,shapeparams,visual)

% If a visualization function has been specified in the visual property of
% the system, use this to draw the body of the system. Otherwise, choose a
% default drawing program based on the system geometry
if exist('visual','var') && isstruct(visual) && isfield(visual,'drawing_function')
    drawing_generator = visual.drawing_function;
else
    % Identify what kind of system is being calculated, and use this to specify how
    % the local connection should be generated
    switch geometry.type

        case {'curvature basis','curvature bases','general curvature'}
            drawing_generator = @fat_backbone;

        case {'n-link chain'}
            drawing_generator = @fat_chain;
            
        case {'branched chain'}
            drawing_generator = @fat_branched_chain;

    end
end

% Evaluate the drawing generator
B = drawing_generator(geometry,shapeparams);