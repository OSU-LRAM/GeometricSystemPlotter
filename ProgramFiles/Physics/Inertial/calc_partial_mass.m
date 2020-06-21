function dM_alphadalpha = calc_partial_mass(s,shapelist)
% Inputs:
%   s: sysplotter system struct for an inertial system
%   shapelist: cell containing n values where each value represents the
%       system's shape at which dM_alphadalpha should be calculated
% Outputs:
%   dM_alphadalpha: Partial derivative of mass with respect to shape
%   variables; cell with n matrices, one for each shape variable

    % Preallocate output
    dM_alphadalpha = cell(size(shapelist));
    % Interpolate in the precalculated dM_alphadalpha grid based on values
    % in shapelist
    for i = 1:length(shapelist)
        dM_alphadalpha{i} = cellfun(@(C) interpn(s.grid.coriolis_eval{:},C,...
            shapelist{:},'spline'),s.coriolisfield.coriolis_eval.content.dM_alphadalpha{i});
    end
end