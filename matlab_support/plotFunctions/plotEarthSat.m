function plotEarthSat(UserPosition, SatPosition)
% PLOTEARTHANDSATELLITES Plots the Earth with longitude and latitude lines,
% user position, and satellite positions in 3D.
%
% INPUTS:
%   - UserPosition: (1x3 vector) User position in ECEF coordinates (meters).
%   - SatPosition: (Nx3 matrix) Satellite positions in ECEF coordinates (meters), where N is the number of satellites.
%
% OUTPUT:
%   A 3D plot of the Earth with longitude and latitude lines, user position, and satellite positions.

    %% Earth Parameters
    earthRadius = 6371e3; % Earth's radius in meters

    %% Create Earth as a sphere
    [xEarth, yEarth, zEarth] = sphere(100); % Generate points for a sphere
    xEarth = xEarth * earthRadius;
    yEarth = yEarth * earthRadius;
    zEarth = zEarth * earthRadius;

    %% Plot the Earth
    figure;
    hold on;
    surf(xEarth, yEarth, zEarth, 'EdgeColor', 'none'); % Plot Earth
    colormap('winter'); % Set Earth color
    alpha(0.7); % Make the Earth slightly transparent

    %% Add Longitude Lines
    lon = -180:30:180; % Longitudes
    lat = -90:1:90; % Latitudes for the lines
    for i = 1:length(lon)
        [x, y, z] = sph2cart(deg2rad(lon(i)), deg2rad(lat), earthRadius);
        plot3(x, y, z, 'k', 'LineWidth', 0.5);
    end

    %% Add Latitude Lines
    lat = -90:30:90; % Latitudes
    lon = -180:1:180; % Longitudes for the lines
    for i = 1:length(lat)
        [x, y, z] = sph2cart(deg2rad(lon), deg2rad(lat(i)), earthRadius);
        plot3(x, y, z*ones(1, length(lon)), 'k', 'LineWidth', 0.5);
    end

    %% Plot User Position
    plot3(UserPosition(1), UserPosition(2), UserPosition(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
    text(UserPosition(1), UserPosition(2), UserPosition(3), ' User', 'VerticalAlignment', 'bottom', 'FontSize', 10);

    %% Plot Satellite Positions
    plot3(SatPosition(:, 1), SatPosition(:, 2), SatPosition(:, 3), 'b^', 'MarkerSize', 8, 'LineWidth', 1.5);
    for i = 1:size(SatPosition, 1)
        text(SatPosition(i, 1), SatPosition(i, 2), SatPosition(i, 3), sprintf(' Sat %d', i), ...
            'VerticalAlignment', 'bottom', 'FontSize', 10);
    end

    %% Customize Plot
    axis equal;
    grid on;
    xlabel('X (meters)');
    ylabel('Y (meters)');
    zlabel('Z (meters)');
    title('User and Satellite Positions with Earth in ECEF Coordinates');
    legend('Earth', 'Longitude & Latitude Lines', 'User Position', 'Satellite Positions', 'Location', 'best');
    view(45, 30); % Set a good viewing angle
end
