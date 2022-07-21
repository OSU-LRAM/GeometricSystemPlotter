function [lagr,lambda,dlambda] = evaluate_lagrangian(objct,const)
%%%%%%%%%%%%%
    % This function calculates the lagrangian function and its derivative
    % for step-optimizer. Some inputs are data structures containing the
    % gradient(gf) and Hessian(hf).
    %
    % Inputs:
    %
    % objct : Data structure containing the gradient and Hessian of
    %   Objective function.
    % const : Data structure containing the gradient and Hessian of
    %   Constraint function.
    %
    % Outputs:
    %
    % lagr : Data structure containing the gradient and Hessian of 
    %   the lagrangian function 
    % lambda : Lagrange multiplier is derived by the moore-penrose
    %   inverse of constraint function.
    % dlambda : The derivative of Lagrange multiplier with respect to the
    %   gait parameter.
    %%%%%%%%%%%%%

    lagr = struct();
    
    % Calculate the lagrange multiplier at the current fourier coefficient.
    const_gf_pinv = pinv(const.gf);
    lambda = const_gf_pinv*objct.gf;
    if size(const.gf,2) == 1
        dlambda = (1/norm(const.gf).^2)*(objct.hf*const.gf+const.hf*objct.gf)...
            -(1/norm(const.gf).^4)*(2*const.hf*const.gf*(const.gf.'*objct.gf));
        
    elseif size(const.gf,2) >= 2
        dlambda = const_gf_pinv*const_gf_pinv.'*const.hf*...
            (eye(length(const.gf))-const.gf*const_gf_pinv)*objct.gf...
            +(eye(length(const.gf))-const.gf*const_gf_pinv)*const.hf...
            *const_gf_pinv.'*const_gf_pinv*objct.gf...
            -const_gf_pinv*const.hf*const_gf_pinv*objct.gf...
            +const_gf_pinv*(objct.hf);
    else
        error("There might be no constraint function.")
    end

    lagr.gf = objct.gf - lambda*const.gf;
    lagr.hf = objct.hf-lambda*const.hf-const.gf*dlambda.';

end