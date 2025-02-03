function R_MCMF_to_ECI = getRmcmf2eci(moonPosECI)
    % computeRMatrix - Computes the rotation matrix from MCMF to ECI
    %
    % Inputs:
    %   moonPosECI - 3x1 vector of the Moon's position in the ECI frame
    %
    % Outputs:
    %   R_MCMF_to_ECI - 3x3 rotation matrix from MCMF to ECI

    % Validate input
    if ~isvector(moonPosECI) || length(moonPosECI) ~= 3
        error('moonPosECI must be a 3x1 vector.');
    end

    % Normalize the Moon's position vector for the MCMF X-axis
    MCMF_X_ECI = -moonPosECI / norm(moonPosECI);

    % Define the Moon's rotation axis inclination (in radians)
    delta = deg2rad(-21.9); % Convert inclination to radians

    % Define the MCMF Z-axis
    MCMF_Z_ECI = [0; sin(delta); cos(delta)];

    % Compute the MCMF Y-axis using the right-hand rule
    MCMF_Y_ECI = cross(MCMF_Z_ECI, MCMF_X_ECI);
    MCMF_Y_ECI = MCMF_Y_ECI / norm(MCMF_Y_ECI);

    % Recompute the MCMF Z-axis to ensure orthogonality
    MCMF_Z_ECI = cross(MCMF_X_ECI, MCMF_Y_ECI);

    % Assemble the rotation matrix
    R_MCMF_to_ECI = [MCMF_X_ECI, MCMF_Y_ECI, MCMF_Z_ECI];
end
