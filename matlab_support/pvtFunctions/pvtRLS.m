function estRxPVT = pvtRLS(refPos, satPos, fracRange, settings, constants)
    %% Load configuration
    num2stepsIterations = settings.rlsIter;
    c = constants.c;
    numSV = size(satPos, 1);

    %% Initialization
    estRxPVT = refPos;
    estRange=zeros(1,numSV);
    %% Iterative estimation
    for kIterations = 1:num2stepsIterations
        % Initialize H matrix
        H = zeros(numSV, 3);

        % Compute range estimates and H matrix
        for kSV = 1:numSV
            estRange(kSV) =  norm(satPos(kSV,:) - estRxPVT(1:3));
            numH = satPos(kSV, :) - estRxPVT(1:3);
            denH = norm(numH);
            H(kSV, 1:3) = -numH / denH;
        end

        % Compute corrections
        estFracRange = rem(estRange / c, 1e-3) * c;
        corrP = fracRange - estFracRange';

        % Update position using least squares
        deltaPVT = ((H' * H) \ H') * corrP;
        estRxPVT = estRxPVT + deltaPVT.';
    end

    % Debug: Calculate position error (optional)
    % PosErrLS = norm(estRxPVT(1:3) - UserPosition); % Can print or store this if needed
end
