function posECI = moonLlh2eci(userLlh, moonPos)
    % userLlhToEci - Converts user geodetic coordinates (LLH) to ECI position
    %
    % Inputs:
    %   userLlh  - 3x1 vector of user latitude (rad), longitude (rad), height (m)
    %   moonPos  - Nx3 matrix of Moon's position in ECI for each update (N rows)
    %   iUpdate  - Index of the current time update (1-based index)
    %
    % Outputs:
    %   posECI   - 3x1 vector of position in the ECI frame

    % Step 1: Convert from LLH to MCMF position
    posMCMF = moonLlh2mcmf(userLlh(1), userLlh(2), userLlh(3));

    % Step 2: Get the rotation matrix from MCMF to ECI
    R_MCMF_to_ECI = getRmcmf2eci(moonPos);

    % Step 3: Convert from MCMF to ECI position
    posECI = mcmf2eci(posMCMF, moonPos, R_MCMF_to_ECI);
end

