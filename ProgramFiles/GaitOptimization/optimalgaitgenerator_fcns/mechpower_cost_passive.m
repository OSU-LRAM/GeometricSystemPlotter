function dcost = mechpower_cost_passive(M,dM,shape,dshape,ddshape,metric)
% Calculates the incremental cost for an inertial system where cost is mechanical power.
% Inputs:
%   M: Mass matrix
%   dM_alphadalpha: Partial of mass matrix with respect to shape variables;
%       must be cell where the ith entry is the partial with respect to the
%       ith shape variable
%   shape: Current shape configuration of system, at which M and
%       dM_alphadalpha were evaluated
%   dshape: Current shape velocity of system
%   ddshape: Current shape acceleration of system

    % Start by calculating the coriolis matrix
    C = calc_coriolis_matrix(dM,shape,dshape);
    % Calculate the torque for this instant of time
    dtau = M*ddshape(:) + C;
    dcost = abs(dtau(2)*dshape(2));
end