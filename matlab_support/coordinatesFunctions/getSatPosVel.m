function [satPos, satVel] = getSatPosVel(constellation, settings, gpsTime, eph)
    % getSatPosVel - Computes satellite positions and velocities for a GNSS constellation
    %
    % Inputs:
    %   constellation - GNSS constellation ('G' for GPS, 'R' for GLONASS, 'E' for Galileo)
    %   settings      - Structure containing settings like numOfUpdates, numOfSatGps, etc.
    %   gpsTime       - Array of GPS time epochs for each update
    %   eph           - Ephemeris data for the selected GNSS constellation
    %
    % Outputs:
    %   satPos - Matrix of satellite positions (Nx3)
    %   satVel - Matrix of satellite velocities (Nx3)

    % Set parameters based on constellation
    switch constellation
        case 'G'
            numberOfSat = settings.numOfSatGps;
            totalNumOfSat = 32;
            prn = settings.prnGps;
        case 'R'
            numberOfSat = settings.numOfSatGlo;
            totalNumOfSat = 24;
            prn = settings.prnGlo;
        case 'E'
            numberOfSat = settings.numOfSatGal;
            totalNumOfSat = 36;
            prn = settings.prnGal;
    end

    % Initialize matrices for satellite positions and velocities
    satPos = zeros(settings.numOfUpdates * numberOfSat, 3);
    satVel = zeros(settings.numOfUpdates * numberOfSat, 3);

    % Initialize waitbar
    hWaitBar = waitbar(0, 'Computing satellite positions and velocities...');

    % Loop through updates to compute positions and velocities
    for iUpdate = 1:settings.numOfUpdates
        % Update the waitbar
        waitbar(iUpdate / settings.numOfUpdates, hWaitBar);

        % Determine the interval for the current update
        intervalGps = (iUpdate - 1) * numberOfSat + 1:iUpdate * numberOfSat;

        % Call satellite_positions to get positions and velocities
        [~, ~, xs, vs, ~, ~, ~] = satellite_positions( ...
            gpsTime(iUpdate), ones(totalNumOfSat, 36), (1:totalNumOfSat), ...
            eph, [], [], zeros(totalNumOfSat, 1), zeros(totalNumOfSat, 1), 0);

        % Assign results to the appropriate intervals
        satPos(intervalGps, :) = xs(prn, :);
        satVel(intervalGps, :) = vs(prn, :);
    end

    % Close the waitbar
    close(hWaitBar);
end
