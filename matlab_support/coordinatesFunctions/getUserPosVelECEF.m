function [userPosECEF, userVelECEF] = getUserPosVelECEF(userLlh, moonPos, utcTime0, settings)
    % getUserPosVelECEF - Computes user positions and velocities in ECEF frametest
    %
    % Inputs:
    %   userLlh      - 3x1 vector of user latitude (rad), longitude (rad), and height (m)
    %   moonPos      - Nx3 matrix of Moon's position in ECI for each update (N rows)
    %   utcTime0     - UTC time of the first update
    %   settings     - Structure containing settings (e.g., numOfUpdates)
    %
    % Outputs:
    %   userPosECEF  - Nx3 matrix of user positions in ECEF frame for all updates
    %   userVelECEF  - Nx3 matrix of user velocities in ECEF frame for all updates

    % Validate inputs
    if size(moonPos, 2) ~= 3
        error('moonPos must be an Nx3 matrix of Moon positions in ECI.');
    end

    if settings.numOfUpdates > size(moonPos, 1)
        error('The number of updates exceeds the number of Moon positions provided.');
    end

    % Initialize matrices for position and velocity
    userPosECEF = zeros(settings.numOfUpdates, 3);
    userVelECEF = zeros(settings.numOfUpdates, 3);

    % Initialize the waitbar
    hWaitBar = waitbar(0, 'Computing user positions and velocities in ECEF frame...');

    % Loop through updates to compute ECEF positions
    for iUpdate = 1:settings.numOfUpdates
        % Convert user LLH to ECI for the current Moon position
        userPosECI = moonLlh2eci(userLlh, moonPos(iUpdate, :)');

        % Convert ECI to ECEF using the given UTC time
        userPosECEF(iUpdate, :) = eci2ecef(utcTime0, userPosECI)';

        % Update the waitbar
        waitbar(iUpdate / settings.numOfUpdates, hWaitBar);
    end

    % Estimate velocities using finite differences
    for iUpdate = 2:settings.numOfUpdates - 1
        % Central difference formula
        userVelECEF(iUpdate, :) = (userPosECEF(iUpdate + 1, :) - userPosECEF(iUpdate - 1, :)) / (2 * settings.moonRate);
    end

    % Forward and backward differences for the first and last updates
    userVelECEF(1, :) = (userPosECEF(2, :) - userPosECEF(1, :)) / settings.moonRate;
    userVelECEF(end, :) = (userPosECEF(end, :) - userPosECEF(end - 1, :)) / settings.moonRate;

    % Close the waitbar
    close(hWaitBar);
end
