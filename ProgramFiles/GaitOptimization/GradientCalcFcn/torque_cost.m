function dcost = torque_cost(M,dM,shape,dshape,ddshape,metric)
% Calculates the incremental cost for an inertial system where cost is torque squared.
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
    % Calculate the torque for this instant of time and return the inner
    % product of the torque with itself
    dtau = M*ddshape(:) + C;
    dcost = dtau'*dtau;
end