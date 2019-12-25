function valid = verify_shape_equations(p,err_tol)
% Verifies if the definition of the shape equations phi_def, dphi_def, and
% ddphi_def are valid against each other by checking the numerical
% derivatives of each against the analytic derivative

% Set a default error tolerance if one isn't provided
if ~exist('err_tol','var') || err_tol == []
    err_tol = 1e-3;
end

% Set the time delta to be used to check the derivatives
dt = 0.001;
% Test ten different values
for i = 1:10
    % Pick a random time value in the interval (0,pi)
    t = pi*rand(1);
    % Take the shape value at the time t and t+dt
    shape = [p.phi_def(t); p.phi_def(t+dt)];
    dshape = [p.dphi_def(t); p.dphi_def(t+dt)];
    ddshape = p.ddphi_def(t);

    err_dshape = dshape(1,:) - (shape(2,:)-shape(1,:))/dt
    err_ddshape = ddshape - (dshape(2,:)-dshape(1,:))/dt

    if any(abs(err_dshape) > err_tol) || any(abs(err_ddshape) > err_tol)
        valid = false;
        return;
    end
end
valid = true;