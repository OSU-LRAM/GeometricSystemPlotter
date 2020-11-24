function ddM_alphadalpha = calc_second_partial_mass(s,shapelist)
% Inputs:
%   s: sysplotter system struct for an inertial system
%   shapelist: cell containing n values where each value represents the
%       system's shape at which ddM_alphadalpha should be calculated
% Outputs:
%   ddM_alphadalpha: Second partial derivative of mass with respect to shape
%   variables; n-by-n cell of n-by-n matrices, one for each shape variable
    
    % Preallocate output
    ddM_alphadalpha = cell(length(shapelist));
    % Interpolate in the precalculated dM_alphadalpha grid based on values
    % in shapelist
    for i = 1:numel(ddM_alphadalpha)
        ddM_alphadalpha{i} = cellfun(@(C) interpn(s.grid.coriolis_eval{:},C,...
            shapelist{:},'spline'),s.coriolisfield.coriolis_eval.content.ddM_alphadalpha{i});
    end
end