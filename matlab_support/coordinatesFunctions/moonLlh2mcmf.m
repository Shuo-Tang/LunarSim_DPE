function posMCMF = moonLlh2mcmf(lat, lon, h)
    % LLHTOMCMF Converts lunar LLH coordinates to MCMF Cartesian coordinates.
    % Inputs:
    %   longitude - Longitude in degrees
    %   latitude  - Latitude in degrees
    %   height    - Height above lunar surface in meters
    % Output:
    %   mcmfCoords - 3x1 vector [X; Y; Z] in MCMF coordinates (meters)

    % Constants
    moonRadius = 1737e3; % Mean radius of the Moon in meters

    % Convert degrees to radians
    latRad = deg2rad(lat);
    lonRad = deg2rad(lon);

    % Compute Cartesian coordinates
    r = moonRadius + h;
    X = r * cos(lonRad) * cos(latRad);
    Y = r * cos(lonRad) * sin(latRad);
    Z = r * sin(lonRad);
    posMCMF = [X; Y; Z];
end