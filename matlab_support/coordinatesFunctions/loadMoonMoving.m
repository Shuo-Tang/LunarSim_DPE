function [moonPos, moonVel] = loadMoonMoving(fileName)
    % loadMoonMoving - Load CSV file with Moon's position and velocity data
    %
    % Inputs:
    %   fileName - Name or path of the CSV file to load
    %
    % Outputs:
    %   moonData - A structure containing:
    %       position - Nx3 matrix of Moon's position (x, y, z)
    %       velocity - Nx3 matrix of Moon's velocity (vx, vy, vz)

    % Check if the file exists
    if ~isfile(fileName)
        error('File "%s" does not exist.', fileName);
    end

    % Load the CSV data
    rawData = readmatrix(fileName);

    % Validate the format of the data
    if size(rawData, 2) ~= 6
        error('The CSV file must contain exactly 6 columns: [x, y, z, vx, vy, vz].');
    end

    % Extract position and velocity components
    moonPos = rawData(:, 1:3); % Columns 1-3: x, y, z
    moonVel = rawData(:, 4:6); % Columns 4-6: vx, vy, vz
end