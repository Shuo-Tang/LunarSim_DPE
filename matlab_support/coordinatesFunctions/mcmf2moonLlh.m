function userLlh = mcmf2moonLlh(posMCMF)
    % mcmf2moonLlh - Converts MCMF Cartesian coordinates to lunar LLH
    %
    % Inputs:
    %   posMCMF - 3x1 vector [X; Y; Z] in MCMF coordinates (meters)
    %
    % Outputs:
    %   userLlh - 3x1 vector of latitude (rad), longitude (rad), height (m)

    % Constants
    moonRadius = 1737e3; % Mean radius of the Moon in meters

    % Compute longitude and latitude
    lat = atan2(posMCMF(2), posMCMF(1));
    lon = atan2(posMCMF(3), sqrt(posMCMF(1)^2 + posMCMF(2)^2));

    % Compute height above lunar surface
    r = norm(posMCMF);
    h = r - moonRadius;

    % Store results
    userLlh = [lat/pi*180; lon/pi*180; h];
end