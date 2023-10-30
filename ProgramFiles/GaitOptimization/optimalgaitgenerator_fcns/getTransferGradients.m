%Gets a cell array of gradients of passive/active fourier coefficients with
%respect to the active coefficients. For N active coefficients (including
%a0 and w), this produces an Nx1 vector of cells.  Each entry i in the cell
%vector is an Nx2 matrix with the gradients of the passive/active coeffs
%wrt the ith active fourier coefficient.  Calculated using the emergent
%transfer function given the current active and passive fourier
%coefficients.
function transferGrads = getTransferGradients(activeCoeffList,passiveCoeffList)

    %Number of fourier coefficients for the active joint
    nCoeffs = numel(activeCoeffList);
    %Number of fourier modes fit (e.g. 4 for a fourier4 fit)
    nModes = (size(activeCoeffList,1)-2)/2;
    %Frequency of the active joint
    w = activeCoeffList(end);
    %Used to avoid division by zero
    small = 1e-10;

    %Reshape coeffs to be an Mx2 matrix for M modes
    activeCoeffs = reshape(activeCoeffList(2:end-1),[2,nModes])';
    passiveCoeffs = reshape(passiveCoeffList(2:end-1),[2,nModes])';

    %Get mode amplitudes and phases from fourier coefficients
    activeMags = zeros(nModes,1);
    passiveMags = zeros(nModes,1);
    activePhases = zeros(nModes,1);
    passivePhases = zeros(nModes,1);
    for i = 1:nModes
        activeMags(i) = norm(activeCoeffs(i,:));
        passiveMags(i) = norm(passiveCoeffs(i,:));
        activePhases(i) = atan2(activeCoeffs(i,2),activeCoeffs(i,1));
        passivePhases(i) = atan2(passiveCoeffs(i,2),passiveCoeffs(i,1));
    end

    %Estimate Output/Input amplitude ratios and Phase Shifts (transfer
    %function values)
    maxRatio = 2;
    magRatios = zeros(nModes,1);
    phaseShifts = zeros(nModes,1);
    for i = 1:nModes
        magRatios(i) = min(passiveMags(i)/max(small,activeMags(i)),maxRatio);
        phaseShifts(i) = passivePhases(i)-activePhases(i);
    end

    %Calculate amplitude-weighted values of transfer function, since at low
    %modal activations the passive joint response is effectively noise 
    %Use these amplitude-weighted values instead of the actual transfer
    %function estimation when the amplitude is below a certain cutoff
    %threshold (1 angular degree of activity).  Actual used value is a
    %sliding scale from 0 activity to threshold activity, with 0 activity
    %using the amplitude-weighted value and threshold activity using the
    %actual transfer function estimate.  This is to make gradients less
    %jumpy around the threshold
    ampWeightedRatio = sum(activeMags.*magRatios)/sum(activeMags);
    ampWeightedPhaseChange = sum(activeMags.*phaseShifts)/sum(activeMags);
    cutoff = pi/180;
    belowCutoff = zeros(nModes,1);
    effectiveMagRatios = magRatios;
    effectivePhaseShifts = phaseShifts;
    for i = 1:nModes
        if activeMags(i) < cutoff
            belowCutoff(i) = 1;
            tradeoff = activeMags(i)/cutoff;
            effectiveMagRatios(i) = tradeoff*magRatios(i) + (1-tradeoff)*ampWeightedRatio;
            effectivePhaseShifts(i) = tradeoff*phaseShifts(i) + (1-tradeoff)*ampWeightedPhaseChange;
        end
    end

    %Derivative of modal magnitude wrt input coefficient
    dMdAB = zeros(nModes,2);
    %Derivative of phase wrt input coefficient
    dPhasedAB = zeros(nModes,2); 
    %Derivative of output/input ratio wrt input coeff (0 if not
    %amplitude-weighted)
    dRwdAB = zeros(nModes,2);
    %Derivative of phase shift wrt input coeff (0 if not
    %amplitude-weighted)
    dPwdAB = zeros(nModes,2);
    %Useful values for gradient calculation
    Msum = sum(activeMags);
    MRsum = sum(activeMags.*magRatios);
    MPsum = sum(activeMags.*phaseShifts);
    for i = 1:nModes
        dMdAB(i,:) = [activeCoeffs(i,1)/max(activeMags(i),small),activeCoeffs(i,2)/max(activeMags(i),small)];
        dPhasedAB(i,:) = [-activeCoeffs(i,2)/(max(activeMags(i),small)^2),activeCoeffs(i,1)/(max(activeMags(i),small)^2)];
        if belowCutoff(i)
            if activeMags(i) < small
                dMdAB(i,:) = [1,1];
            end
            dRavgdM = (magRatios(i)*Msum - MRsum)/(Msum^2);
            dPavgdM = (phaseShifts(i)*Msum - MPsum)/(Msum^2);
            dRwdM = magRatios(i)/cutoff - ampWeightedRatio/cutoff + (1-magRatios(i)/cutoff)*dRavgdM;
            dPwdM = phaseShifts(i)/cutoff - ampWeightedPhaseChange/cutoff + (1-magRatios(i)/cutoff)*dPavgdM;
            dRwdAB(i,:) = dRwdM*dMdAB(i,:);
            dPwdAB(i,:) = dPwdM*dMdAB(i,:);
        end
    end

    %Calculate gradients of passive coefficients with respect fo active
    %ones
    %Derivative of first passive coefficient in mode wrt both passive and
    %active coefficient in active mode
    dAoutdAB = zeros(nModes,2);
    %Derivative of second passive coefficient in mode wrt both passive and
    %active coefficient in active mode
    dBoutdAB = zeros(nModes,2);
    for i = 1:nModes
        dAoutdAB(i,:) = dRwdAB(i,:)*activeMags(i)*cos(activePhases(i)+effectivePhaseShifts(i)) - ...
            dPwdAB(i,:)*effectiveMagRatios(i)*activeMags(i)*sin(activePhases(i)+effectivePhaseShifts(i)) + ...
            dMdAB(i,:)*effectiveMagRatios(i)*cos(activePhases(i)+effectivePhaseShifts(i)) - ...
            dPhasedAB(i,:)*effectiveMagRatios(i)*activeMags(i)*sin(activePhases(i)+effectivePhaseShifts(i));

        dBoutdAB(i,:) = dRwdAB(i,:)*activeMags(i)*sin(activePhases(i)+effectivePhaseShifts(i)) + ...
            dPwdAB(i,:)*effectiveMagRatios(i)*activeMags(i)*cos(activePhases(i)+effectivePhaseShifts(i)) + ...
            dMdAB(i,:)*effectiveMagRatios(i)*sin(activePhases(i)+effectivePhaseShifts(i)) + ...
            dPhasedAB(i,:)*effectiveMagRatios(i)*activeMags(i)*cos(activePhases(i)+effectivePhaseShifts(i));
    end

    %Repackage gradients in a convenient form for manipulation in optimizer
    transferGrads = cell(nCoeffs,1);
    for i = 1:nCoeffs
        transferGrads{i} = zeros(nCoeffs,2);
        transferGrads{i}(i,2) = 1;
        if i == nCoeffs
            transferGrads{i}(i,1) = 1;
        elseif i ~= 1
            modeNum = floor(i/2);
            Aind = 2*modeNum;
            Bind = Aind+1;
            if ~mod(i,2)
                transferGrads{i}(Aind,1) = dAoutdAB(modeNum,1);
                transferGrads{i}(Bind,1) = dBoutdAB(modeNum,1);
            else
                transferGrads{i}(Aind,1) = dAoutdAB(modeNum,2);
                transferGrads{i}(Bind,1) = dBoutdAB(modeNum,2);
            end
        end
    end

end

