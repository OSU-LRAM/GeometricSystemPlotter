function pdot=ode_step_optimal(p,s,n,dimension,direction,lb,ub,nfparam,method)%,lb,ub,writerObj)
    %%%%%%%%%%%%%
    % This function calculates efficiency (or displacement, if
    % that is the objective function) and its gradient with respect to the coefficients obtained
    % by the fourier series parametrization
    %
    % Inputs:
    %
    % p: Matrix containing the Fourier series coefficients
    % s: System file which contains the connection vector field, CCF's and
    %   metric data
    % dimension: Indicates the number of shape variables of the system.
    % n: The number of points desired in a direct transcription parametrization
    %   of the gaits
    % direction: direction in which to optimize motion: 1-x, 2-y, 3-theta,
    %   4-steering gait(x and theta)
    % costfunction: costfunction type to optimize for
    % lb: Lower bound of shape variables for each point which is obtained from the grid inside
    %   which an optimal gait is desired
    % ub: Upper bound of shape variables for each point which ll'l;is obtained from the grid inside
    %   which an optimal gait is desired
    % method: 'Hessian' - Based on the null vector of Hessian, the result
    %   vector is evaluated. Since the vector always point to the next step
    %   gait(stable), ode45 tend to skip some steps.
    %   'Perturb' - Based on the stable point of the gradient of
    %   Lagrangian function, the result vector is evaluated. At first, the
    %   vector does not point to the next step gait but it will fall into
    %   the stable point(next step) eventually. So, ode45 will try to
    %   calculate every step's solution.
    %
    % Outputs:
    %
    % pdot : the resulting vector is the negative gradient of displacement
    %   projected onto the null space of Hessian.
    %%%%%%%%%%%%%

    % Because of ODE solver, size of y is the number of fourier coefficient*dimension x 1
    p=reshape(p,[nfparam dimension]);

    disp = struct();
    stroke = struct();
    eff = struct();

    % when not optimizing the steering gait, consider only one direction's
    % Hessian and gradient
    if direction ~= 4
        [jacobfourier,temp_disp,temp_stroke] = evaluate_jacobian_fourier(p,s,n,dimension,direction,lb,ub);
        disp.f = temp_disp;
        stroke.f = temp_stroke;
        if strcmpi(method,'Hessian')
            % TODO : evaluate the hessian mathmatically.
            hessfourier = numelhessfourier(p,s,n,dimension,direction,lb,ub,jacobfourier);
            disp.hf = cell2mat(hessfourier.disp);
            stroke.hf = cell2mat(hessfourier.stroke);
            repuls.hf = cell2mat(hessfourier.repuls);
        end
        % when optimizing the steering gait, consider x and theta direction's
        % Hessian and gradient
    else
        [jacobfourierx,temp_disp,temp_stroke] = evaluate_jacobian_fourier(p,s,n,dimension,1,lb,ub);
        disp.f = temp_disp;
        stroke.f = temp_stroke;
        [jacobfourier,~,~] = evaluate_jacobian_fourier(p,s,n,dimension,3,lb,ub);
        % Calculate the Hessian of efficiency in the x direction.
        if strcmpi(method,'Hessian')
            hessfourierx = numelhessfourier(p,s,n,dimension,1,lb,ub,jacobfourierx);
            hessfourier = numelhessfourier(p,s,n,dimension,3,lb,ub,jacobfourier);
            disp.hf = cell2mat(hessfourierx.disp);
            strokex.hf = cell2mat(hessfourierx.stroke);
            disp.hf = cell2mat(hessfourier.disp);
            stroke.hf = cell2mat(hessfourier.stroke);
        end
    end
    %% Record the fourier coefficient at step-optimal gaits.
    global Eff currentDisp bestDisp currentCost stepOptimalGaits;
    if direction == 4
        global rotOptMod
        dir = 1;
    else
        dir = direction;
    end

    % Calculate the ratio the current displacement to the maximum
    % displacement(the seed gait). If the condition at each step is
    % satisfied, store the gait information to stepOptimalGaits.
    curDispPer= disp.f(dir)/bestDisp;
    prevDispPer = currentDisp(dir)/bestDisp;
    currentCost = stroke.f;

    CheckPoint = [3/4, 2/4, 1/4];
    for i = 1:3
        if (abs(curDispPer - CheckPoint(i)) < abs(prevDispPer - CheckPoint(i)))...
                && (abs(curDispPer - CheckPoint(i)) < 1e-2)
            stepOptimalGaits{i+1,1} = reshape(p,[nfparam dimension]);
            stepOptimalGaits{i+1,2} = [disp.f(dir),stroke.f];
            stepOptimalGaits{i+1,3} = [jacobfourier.disp jacobfourier.stroke];
        end
    end

    currentDisp = disp.f;
    Eff = [Eff; currentDisp(3)/currentCost currentDisp(1)/currentCost];

    %% calcuating pdot so that lagrange equation is zero.

    % Reshape jacobdisp and jacobstorke for ODE Solver.
    if direction ~= 4
        disp.gf = reshape(jacobfourier.disp,[(nfparam-1)*dimension 1]);
        stroke.gf = reshape(jacobfourier.stroke,[(nfparam-1)*dimension 1]);
        repuls.gf = reshape(jacobfourier.repuls,[(nfparam-1)*dimension 1]);
        eff.gf = reshape(jacobfourier.eff,[(nfparam-1)*dimension 1]);        
    else
        disp.gf = reshape(jacobfourier.disp,[(nfparam-1)*dimension 1]);
        stroke.gf = reshape(jacobfourier.stroke,[(nfparam-1)*dimension 1]);
        eff.gf  = reshape(jacobfourier.eff,[(nfparam-1)*dimension 1]);

        dispx.gf  = reshape(jacobfourierx.disp,[(nfparam-1)*dimension 1]);
        strokex.gf = reshape(jacobfourierx.stroke,[(nfparam-1)*dimension 1]);
        effx.gf = reshape(jacobfourierx.eff,[(nfparam-1)*dimension 1]);
    end

    if strcmpi(method,'Hessian')
        if (direction ~= 4)
            lagr = evaluate_lagrangian(stroke,disp);
            optdir = -1;
        else
            eff.hf = fraction_hessian(disp(3),cost,...
                gradfourier.disp,gradfourier.stroke,...
                hessfourier.disp,hessfourier.stroke);
            effx.hf = fraction_hessian(disp(1),cost,...
                gradfourierx.disp,gradfourierx.stroke,...
                hessfourierx.disp,hessfourierx.stroke);

            if (rotOptMod == 0)
                lagr = evaluate_lagrangian(effx,eff);
                optdir = -1;
            else
                lagr = evaluate_lagrangian(eff,effx);
                optdir = 1;
            end
        end

        % The right null space of matrix is the columns of V corresponding
        % to singular values equal to zero.
        [~,S,V]=svd(lagr.hf);
        S=diag(real(S));

        pdot = zeros(length(S),1);

        % Find the null space of Hessian. Its tolerance is 0.1.
        Sidx = find((S < 0.1) & (S > -0.1)).';
        
        % Project the gradient vector onto the null space
        for i = Sidx
            pdot=pdot+optdir*(V(:,i).'*disp.gf)/norm(V(:,i))^2*V(:,i);
        end
        if isnan(pdot)
            pdot = zeros(size(pdot));
        end

    elseif strcmpi(method,'Perturb')
        % The solution vector is the gradient of the Lagrangian equation
        % with the offset of the lagrange multiplier.
        pdot = -stroke.gf+(1-0.02)*pinv(disp.gf)...
            *stroke.gf*disp.gf;
    end

    % The frequency should not change.
    pdot = reshape(pdot,[nfparam - 1 dimension]);
    pdot = [pdot; zeros(1,dimension)];
    pdot = reshape(pdot,[nfparam*dimension 1]);


end