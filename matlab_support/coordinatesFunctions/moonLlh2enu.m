function menuCoords = moonLlh2enu(userLlh, refLlh)
    % moonLlh2menu - Converts lunar LLH coordinates to local menu coordinates
    %
    % Inputs:
    %   userLlh - 3x1 vector of user latitude (rad), longitude (rad), height (m)
    %   refLlh  - 3x1 vector of reference latitude, longitude, height (m)
    %
    % Outputs:
    %   menuCoords - 3x1 vector of position in local menu coordinates (meters)

    % Convert reference LLH to MCMF
    refMCMF = moonLlh2mcmf(refLlh(1), refLlh(2), refLlh(3));
    
    % Convert user LLH to MCMF
    userMCMF = moonLlh2mcmf(userLlh(1), userLlh(2), userLlh(3));
    
    % Compute the relative position in MCMF frame
    deltaMCMF = userMCMF - refMCMF;
    
    % Compute the local ENU rotation matrix
    lat = refLlh(1);
    lon = refLlh(2);
    R_enu = [-sin(lon), cos(lon), 0;
             -sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat);
              cos(lat)*cos(lon), cos(lat)*sin(lon), sin(lat)];
    
    % Convert to local menu coordinates
    menuCoords = R_enu * deltaMCMF;
end