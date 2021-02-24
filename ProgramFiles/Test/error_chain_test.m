% test script for error_chain.m
% may also be used a reference for use

% test error handling
disp('Errors:')
% no geometry
try
    error_chain();
catch e
    disp(e.message);
end
% empty geometry
geometry = struct;
try
    error_chain(geometry);
catch e
    disp(e.message);
end
% incorrect geometry elements
geometry.elements = {1 2 3 4};
geometry.types = cell(4);
try
    error_chain(geometry);
catch e
    disp(e.message);
end
% incorrect geometry types
geometry.elements = {GroupElement.LINEAR};
geometry.types = {1 2 3 4};
try
    error_chain(geometry);
catch e
    disp(e.message);
end

% construct planar error chain; get 1st order error
geometry.elements = {GroupElement.ROTARY, GroupElement.ROTARY,...
                     GroupElement.LINEAR};
geometry.types = {'j', 'e', 'l'};
configuration.errors = sym('theta');
disp('Expect warning about order:');
err_lae = error_chain(geometry, configuration);
disp('End effector error LAE:');
err_lae.vector



