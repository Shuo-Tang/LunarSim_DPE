function userLlh = eci2moonLlh(posECI, moonPos)
    % eci2moonLlh - Converts ECI position to lunar geodetic coordinates (LLH)
    %
    % Inputs:
    %   posECI   - 3x1 vector of position in the ECI frame
    %   moonPos  - Nx3 matrix of Moon's position in ECI for each update (N rows)
    %
    % Outputs:
    %   userLlh  - 3x1 vector of user latitude (rad), longitude (rad), height (m)

    % Step 1: Get the rotation matrix from ECI to MCMF
    R_ECI_to_MCMF = getRmcmf2eci(moonPos)';  % Inverse of MCMF to ECI

    % Step 2: Convert from ECI to MCMF position
    posMCMF = R_ECI_to_MCMF * (posECI - moonPos);

    % Step 3: Convert from MCMF to LLH
    userLlh = mcmf2moonLlh(posMCMF);
end
