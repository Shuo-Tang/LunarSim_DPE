clear;
clc;
close all;

%% =========== Load Satellite Data =============
satEphTable = readtable("../exportData/satEph.csv");
satPosVelTable = readtable("../exportData/satPosVel.csv");
moonDataPath = "../testData/horizons_results.txt";
% === manually headers according to the system ===
satEphTable.Properties.VariableNames{1} = 'Constellation';
satEphTable.Properties.VariableNames{2} = 'SatId';
satEphTable.Properties.VariableNames{3} = 'GpsTime';
% === convert to an data array ===
satEph = table2array(satEphTable(:, 4:end));
satPosVel = table2array(satPosVelTable(:, 2:7));

%% =========== Constants / Settings =============
mu = 3.986004418e14; % Earth's gravitational constant (m^3/s^2)
earthRadius = 6.356752314245179e6; %6371e3; % Earth's radius (meters)
moonRadius = 1737e3; % Moon's radius (meters)
OMEGAE_DOT_GPS = 7.2921151467e-5;             % GPS     Angular velocity of the Earth rotation [rad/s]
OMEGAE_DOT_GLO = 7.292115e-5;                 % GLONASS Angular velocity of the Earth rotation [rad/s]
OMEGAE_DOT_GAL = 7.2921151467e-5;             % Galileo Angular velocity of the Earth rotation [rad/s]
OMEGAE_DOT_BDS = 7.292115e-5;                 % BeiDou  Angular velocity of the Earth rotation [rad/s]
CIRCLE_RAD = 2* pi;

numOfSat = size(satEph, 1);
t0 = satEphTable.GpsTime(1);  
timeStep = 600;
timeSpan = 20*3600;
T = t0: timeStep: t0 + timeSpan;
utcDate = gps2utcDate(t0);

enablePosLoad = true;
enableEarthPot = true;
enableEarthAxisPlot = true;
enableSatPlot = true;
enableMoonPlot = true;
enableMoonAxisPlot = true;

enableECEFPlot = true;
enableECIPlot = true;

enableTestPlot = false;
enableValidSatPlot = true;
enableGDOPPlot = false;
%% =========== Compute Orbits =============
% === Loop for all satellites ===
satPosECEF = cell(numOfSat, 1);
satPosECI = cell(numOfSat, 1);
satPosCurECEF = zeros(numOfSat, 3);
satPosCurECI = zeros(numOfSat, 3);

for iSat = 1: numOfSat
    constellation = satEphTable.Constellation{iSat};
    Eph = satEph(iSat,:);
    switch constellation
        case 'G'
            Omegae_dot = OMEGAE_DOT_GPS;
        case 'R'
            Omegae_dot = OMEGAE_DOT_GLO;
        case 'E'
            Omegae_dot = OMEGAE_DOT_GAL;
        case 'C'
            Omegae_dot = OMEGAE_DOT_BDS;
    end
    % not GLONASS
    if ~strcmpi(constellation, 'R')
        % Extract orbital elements
        roota = Eph(11);
        ecc = Eph(9);% Eccentricity
        Omega0 = Eph(14);
        i0 = Eph(16); % Inclination (rad)
        Omega = Eph(18); % Argument of perigee (rad)
        Omega_dot = Eph(19);
        M0 = Eph(7); % Mean anomaly (rad)
        cuc = Eph(8);
        cus = Eph(10);
        crc = Eph(17);
        crs = Eph(5);
        cic = Eph(13);
        cis = Eph(15);
        idot = Eph(20);
        toe = Eph(12); % Reference time of ephemeris (s)
        week = Eph(22);
        time_eph = weektow2time(week, toe, constellation);

        %eccentric anomaly
        [Ek, n] = ecc_anomaly_0(T, Eph, constellation);

        cr = 2 * pi;

        A = roota*roota;             %semi-major axis
        tk = check_t(T - time_eph);  %time from the ephemeris reference epoch
        % 
        fk = atan2(sqrt(1-ecc^2)*sin(Ek), cos(Ek) - ecc);    %true anomaly
        phik = fk + Omega;                           %argument of latitude
        phik = rem(phik,cr);

        uk = phik                + cuc*cos(2*phik) + cus*sin(2*phik); %corrected argument of latitude
        rk = A*(1 - ecc*cos(Ek)) + crc*cos(2*phik) + crs*sin(2*phik); %corrected radial distance
        ik = i0 + idot*tk        + cic*cos(2*phik) + cis*sin(2*phik); %corrected inclination of the orbital plane

        % correct longitude of the ascending node
        Omegak = Omega0 ;%+ (Omega_dot - Omegae_dot)*tk - Omegae_dot*toe;
        Omegak = rem(Omegak + CIRCLE_RAD, CIRCLE_RAD);

        % satellite positions in the orbital plane
        x1k = cos(uk).*rk;
        y1k = sin(uk).*rk;

        % Convert to 3D ECEF coordinates
        % For visualization, the ECEF transformation is only available for
        % current time epoch: Omegak = Omega0
        x_ecef = x1k .* cos(Omegak) - y1k .* cos(ik) .* sin(Omegak);
        y_ecef = x1k .* sin(Omegak) + y1k .* cos(ik) .* cos(Omegak);
        z_ecef = y1k .* sin(ik);
        rECI = zeros(3, length(T));
        for iT = 1: length(T)
            rECI(:,iT) = ecef2eci(utcDate, [x_ecef(iT), y_ecef(iT), z_ecef(iT)].');      
        end
        satPosECI{iSat} = rECI / 1e3;

        % Convert to kilometers
        x_ecef = x_ecef / 1e3;
        y_ecef = y_ecef / 1e3;
        z_ecef = z_ecef / 1e3;
        satPosECEF{iSat} = [x_ecef; y_ecef; z_ecef];


        % compute current position
        if enablePosLoad
            x_cur = satPosVel(iSat, 1) / 1e3;
            y_cur = satPosVel(iSat, 2) / 1e3;
            z_cur = satPosVel(iSat, 3) / 1e3;
        else
            % Convert to kilometers
            x_cur = x_ecef(1);
            y_cur = y_ecef(1);
            z_cur = z_ecef(1);
        end
        satPosCurECEF(iSat, :) = [x_cur, y_cur, z_cur]; 
        rCurECI = ecef2eci(utcDate, [x_cur, y_cur, z_cur].' * 1e3);
        satPosCurECI(iSat, :) = rCurECI.' / 1e3;

    % GLONASS
    else
        % GLONASS orbits computation is not provided at this point
        % % Propagate orbits using position and velocity
        % % Extract initial position and velocity
        % initialPos = satPosVel(iSat, 1:3)'; % 3x1 position vector (meters)
        % initialVel = satPosVel(iSat, 4:6)'; % 3x1 velocity vector (meters/second)
        % 
        % % Numerical integration using ode45
        % state = [initialPos; initialVel];
        % options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
        % [~, stateHistory] = ode45(@(t, y) equationsOfMotion(t, y, mu), T, state, options);
        % 
        % % Extract position components and convert to kilometers
        % satPosECEF{iSat} = stateHistory(:, 1:3) / 1e3;
        % % Convert to ECI
        % for iT = 1: length(T)
        %     rECI(:,iT) = ecef2eci(utcDate, [stateHistory(iT, 1), stateHistory(iT, 2), stateHistory(iT, 3)].');      
        % end
        % satPosECI{iSat} = rECI / 1e3;

        % compute current position
        if enablePosLoad
            x_cur = satPosVel(iSat, 1) / 1e3;
            y_cur = satPosVel(iSat, 2) / 1e3;
            z_cur = satPosVel(iSat, 3) / 1e3;
        % else
            % Convert to kilometers
            % x_cur = stateHistory(1, 1) / 1e3;
            % y_cur = stateHistory(1, 2) / 1e3;
            % z_cur = stateHistory(1, 3) / 1e3;
        end
        satPosCurECEF(iSat, :) = [x_cur, y_cur, z_cur]; 
        rCurECI = ecef2eci(utcDate, [x_cur, y_cur, z_cur].' * 1e3);
        satPosCurECI(iSat, :) = rCurECI.' / 1e3;

    end
end


%% =========== Plots =============
if enableECIPlot
    % === Initialize figure ===
    figure(101);
    hold on;
    grid on;
    xlabel('X (km)', 'FontSize', 18);
    ylabel('Y (km)', 'FontSize', 18);
    zlabel('Z (km)', 'FontSize', 18);
    title('Satellite Orbits (ECI)', 'FontSize', 18);
    view(30, 14)
    % === Draw Earth ===
    if enableEarthPot
        [x, y, z] = sphere(50);
        earthP = surf(x * earthRadius / 1e3, y * earthRadius / 1e3, z * earthRadius / 1e3, ...
            'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', [0 0 1], ...
            'DisplayName', 'Earth');
        % Plot Earth's ECI coordinate axes
        x_ECI = [1, 0, 0];
        y_ECI = [0, 1, 0];
        z_ECI = [0, 0, 1];
        scale = 1e5;
        earthXP = quiver3(0, 0, 0, scale * x_ECI(1), scale * x_ECI(2), scale * x_ECI(3), ...
            'r', 'LineWidth', 2, 'DisplayName', 'ECI X-axis');
        earthYP = quiver3(0, 0, 0, scale * y_ECI(1), scale * y_ECI(2), scale * y_ECI(3), ...
            'g', 'LineWidth', 2, 'DisplayName', 'ECI Y-axis');
        earthZP = quiver3(0, 0, 0, scale * z_ECI(1), scale * z_ECI(2), scale * z_ECI(3), ...
            'b', 'LineWidth', 2, 'DisplayName', 'ECI Z-axis');
    end
    % === Draw Satellites Orbits ===
    if enableSatPlot
        for iSat = 1: numOfSat
            % Plot orbit trajectory
            if ~isempty(satPosECI{iSat})
                satOP(iSat) = plot3(satPosECI{iSat}(1,:), satPosECI{iSat}(2,:), satPosECI{iSat}(3,:), ...
                    'DisplayName', satEphTable.SatId{iSat});
            end
            % Plot current position
            satCP(iSat) = plot3(satPosCurECI(iSat, 1), satPosCurECI(iSat, 2), satPosCurECI(iSat, 3), ...
                '.', 'MarkerSize', 20, 'Color', '#0072BD', 'HandleVisibility', 'off');
        end
    end
    % === Draw Moon ===
    if enableMoonPlot
        moonPosECI = getMoonPosECI(utcDate);
        [xm, ym, zm] = sphere(50);
        moonP = surf((xm * moonRadius + moonPosECI(1)) / 1e3 , ...
             (ym * moonRadius + moonPosECI(2)) / 1e3 , ...
             (zm * moonRadius + moonPosECI(3)) / 1e3, ...
             'FaceAlpha', 1, 'EdgeColor', [0.5 0.5 0.5], 'FaceColor', [0.5 0.5 0.5], ...
             'DisplayName', 'Moon');
    end
    if enableMoonAxisPlot
        % Plot Moon's MCMF coordinate axes in ECI
        MCMF_X = - moonPosECI / norm(moonPosECI);
        % MCMF Z-axis: Moon's rotation axis inclined by 21.9° to Earth's equator
        delta = deg2rad(-21.9);
        MCMF_Z = [0; sin(delta); cos(delta)];
        % MCMF Y-axis: Right-hand rule
        MCMF_Y = cross(MCMF_Z, MCMF_X);
        MCMF_Y = MCMF_Y / norm(MCMF_Y);
        MCMF_Z = cross(MCMF_X, MCMF_Y);
    
        % Scale for visualization
        scale = 5e4;
        
        % Plot MCMF axes
        moonXP = quiver3(moonPosECI(1) / 1e3, moonPosECI(2) / 1e3, moonPosECI(3) / 1e3, scale * MCMF_X(1), ...
            scale * MCMF_X(2), scale * MCMF_X(3), 'r', 'LineWidth', 2, ...
            'LineStyle', ':', 'DisplayName', 'MCMF X-axis');
        moonYP = quiver3(moonPosECI(1) / 1e3, moonPosECI(2) / 1e3, moonPosECI(3) / 1e3, scale * MCMF_Y(1), ...
            scale * MCMF_Y(2), scale * MCMF_Y(3), 'g', 'LineWidth', 2, ...
            'LineStyle', ':', 'DisplayName', 'MCMF Y-axis');
        moonZP = quiver3(moonPosECI(1) / 1e3, moonPosECI(2) / 1e3, moonPosECI(3) / 1e3, scale * MCMF_Z(1), ...
            scale * MCMF_Z(2), scale * MCMF_Z(3), 'b', 'LineWidth', 2, ...
            'LineStyle', ':', 'DisplayName', 'MCMF Z-axis');
    end

    % === Draw a test surface point ===
    
    R_MCMF_to_ECI = [MCMF_X, MCMF_Y, MCMF_Z];
    % disp('Rotation Matrix from MCMF to Earth ECI:');
    % disp(R_MCMF_to_ECI);
    lat = 0;
    lon = 90;
    h = 0;
    fprintf("ECI:  The test point is on the Moon at Lat %f, Long %f, h %f\n", lat, lon, h)
    posMCMF = llh2MCMF(lat, lon, h);
    posECI = mcmf2ECI(posMCMF, moonPosECI, R_MCMF_to_ECI);
    if enableTestPlot
        pointP = plot3(posECI(1) / 1e3, posECI(2) / 1e3, posECI(3) / 1e3, '*', 'MarkerSize', 10, ...
                'Color', "#D95319", 'lineWidth', 2, 'DisplayName', 'Test Point');
    end
    % === config the figure ===
    axis equal;
    xlim1 = min([0, moonPosECI(1)]) / 1e3 - 1e5;
    xlim2 = max([0, moonPosECI(1)]) / 1e3 + 1e5;
    ylim1 = min([0, moonPosECI(2)]) / 1e3 - 1e5;
    ylim2 = max([0, moonPosECI(2)]) / 1e3 + 1e5;
    zlim1 = min([0, moonPosECI(3)]) / 1e3 - 1e5;
    zlim2 = max([0, moonPosECI(3)]) / 1e3 + 1e5;
    xlim([xlim1, xlim2])
    ylim([ylim1, ylim2])
    zlim([zlim1, zlim2])
    % legend([earthP, earthXP, earthYP, earthZP, moonP, moonXP, moonYP, moonZP, satOP, pointP]);
    legend([earthP, earthXP, earthYP, earthZP, moonP, moonXP, moonYP, moonZP], 'FontSize', 12);
    hold off;
end

if enableECEFPlot
    % === Initialize figure ===
    figure(102);
    hold on;
    grid on;
    xlabel('X (km)');
    ylabel('Y (km)');
    zlabel('Z (km)');
    title('Satellite Orbits (ECEF)');
    % === Draw Earth ===
    if enableEarthPot
        [x, y, z] = sphere(50);
        earthP = surf(x * earthRadius / 1e3, y * earthRadius / 1e3, z * earthRadius / 1e3, ...
            'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', [0 0 1], ...
            'DisplayName', 'Earth');
        % Plot Earth's ECI coordinate axes
        x1 = [1, 0, 0];
        y1 = [0, 1, 0];
        z1 = [0, 0, 1];
        x_ECI = eci2ecef(utcDate, x1).';
        y_ECI = eci2ecef(utcDate, y1).';
        z_ECI = eci2ecef(utcDate, z1).';

        scale = 1e5;
        earthXP = quiver3(0, 0, 0, scale * x_ECI(1), scale * x_ECI(2), scale * x_ECI(3), ...
            'r', 'LineWidth', 2, 'DisplayName', 'ECI X-axis');
        earthYP = quiver3(0, 0, 0, scale * y_ECI(1), scale * y_ECI(2), scale * y_ECI(3), ...
            'g', 'LineWidth', 2, 'DisplayName', 'ECI Y-axis');
        earthZP = quiver3(0, 0, 0, scale * z_ECI(1), scale * z_ECI(2), scale * z_ECI(3), ...
            'b', 'LineWidth', 2, 'DisplayName', 'ECI Z-axis');
    end
    % === Draw Satellites Orbits ===
    if enableSatPlot
        for iSat = 1: numOfSat
            % Plot orbit trajectory
            if ~isempty(satPosECEF{iSat})
                satOP(iSat) = plot3(satPosECEF{iSat}(1,:), satPosECEF{iSat}(2,:), satPosECEF{iSat}(3,:), ...
                    'DisplayName', satEphTable.SatId{iSat});
            end
            % Plot current position
            satCP(iSat) = plot3(satPosCurECEF(iSat, 1), satPosCurECEF(iSat, 2), satPosCurECEF(iSat, 3), ...
                '.', 'MarkerSize', 20, 'Color', '#0072BD', 'HandleVisibility', 'off');
        end
    end
    % === Draw Moon ===
    if enableMoonPlot
        moonPosECI = getMoonPosECI(utcDate);
        moonPosECEF = eci2ecef(utcDate, moonPosECI);
        [xm, ym, zm] = sphere(50);
        moonP = surf((xm * moonRadius + moonPosECEF(1)) / 1e3 , ...
             (ym * moonRadius + moonPosECEF(2)) / 1e3 , ...
             (zm * moonRadius + moonPosECEF(3)) / 1e3, ...
             'FaceAlpha', 1, 'EdgeColor', [0.5 0.5 0.5], 'FaceColor', [0.5 0.5 0.5], ...
             'DisplayName', 'Moon');
        fprintf("Position of the moon's center in ECI: x = %f, y = %f, z = %f \n", ...
            moonPosECI(1), moonPosECI(2), moonPosECI(3))
        fprintf("Position of the moon's center in ECEF: x = %f, y = %f, z = %f \n", ...
            moonPosECEF(1), moonPosECEF(2), moonPosECEF(3))
    end
    if enableMoonAxisPlot
        % Plot Moon's MCMF coordinate axes in ECEF
        MCMF_X_ECI = - moonPosECI / norm(moonPosECI);
        % MCMF Z-axis: Moon's rotation axis inclined by 21.9° to Earth's equator
        delta = deg2rad(-21.9);
        MCMF_Z_ECI = [0; sin(delta); cos(delta)];
        % MCMF Y-axis: Right-hand rule
        MCMF_Y_ECI = cross(MCMF_Z_ECI, MCMF_X_ECI);
        MCMF_Y_ECI = MCMF_Y_ECI / norm(MCMF_Y_ECI);
        MCMF_Z_ECI = cross(MCMF_X_ECI, MCMF_Y_ECI);
        % Convert to ECEF plots
        MCMF_X = eci2ecef(utcDate, MCMF_X_ECI);
        MCMF_Y = eci2ecef(utcDate, MCMF_Y_ECI);
        MCMF_Z = eci2ecef(utcDate, MCMF_Z_ECI);
        % Scale for visualization
        scale = 5e4;
        
        % Plot MCMF axes
        moonXP = quiver3(moonPosECEF(1) / 1e3, moonPosECEF(2) / 1e3, moonPosECEF(3) / 1e3, scale * MCMF_X(1), ...
            scale * MCMF_X(2), scale * MCMF_X(3), 'r', 'LineWidth', 2, ...
            'LineStyle', ':', 'DisplayName', 'MCMF X-axis');
        moonYP = quiver3(moonPosECEF(1) / 1e3, moonPosECEF(2) / 1e3, moonPosECEF(3) / 1e3, scale * MCMF_Y(1), ...
            scale * MCMF_Y(2), scale * MCMF_Y(3), 'g', 'LineWidth', 2, ...
            'LineStyle', ':', 'DisplayName', 'MCMF Y-axis');
        moonZP = quiver3(moonPosECEF(1) / 1e3, moonPosECEF(2) / 1e3, moonPosECEF(3) / 1e3, scale * MCMF_Z(1), ...
            scale * MCMF_Z(2), scale * MCMF_Z(3), 'b', 'LineWidth', 2, ...
            'LineStyle', ':', 'DisplayName', 'MCMF Z-axis');
    end

    % === Draw a test surface point ===
    
    R_MCMF_to_ECI = [MCMF_X_ECI, MCMF_Y_ECI, MCMF_Z_ECI];
    % disp('Rotation Matrix from MCMF to Earth ECI:');
    % disp(R_MCMF_to_ECI);
    lat = 0;
    lon = 90;
    h = 0;
    fprintf("ECEF: The test point is on the Moon at Lat %f, Long %f, h %f\n", lat, lon, h)
    posMCMF = llh2MCMF(lat, lon, h);
    posECI = mcmf2ECI(posMCMF, moonPosECI, R_MCMF_to_ECI);
    posECEF = eci2ecef(utcDate, posECI);
    if enableTestPlot
        pointP = plot3(posECEF(1) / 1e3, posECEF(2) / 1e3, posECEF(3) / 1e3, '*', 'MarkerSize', 10, ...
                'Color', "#D95319", 'lineWidth', 2, 'DisplayName', 'Test Point');
    end

    % === config the figure ===
    axis equal;
    xlim1 = min([0, moonPosECEF(1)]) / 1e3 - 1e5;
    xlim2 = max([0, moonPosECEF(1)]) / 1e3 + 1e5;
    ylim1 = min([0, moonPosECEF(2)]) / 1e3 - 1e5;
    ylim2 = max([0, moonPosECEF(2)]) / 1e3 + 1e5;
    zlim1 = min([0, moonPosECEF(3)]) / 1e3 - 1e5;
    zlim2 = max([0, moonPosECEF(3)]) / 1e3 + 1e5;
    xlim([xlim1, xlim2])
    ylim([ylim1, ylim2])
    zlim([zlim1, zlim2])
    % legend([earthP, earthXP, earthYP, earthZP, moonP, moonXP, moonYP, moonZP, satOP, pointP]);
    legend([earthP, earthXP, earthYP, earthZP, moonP, moonXP, moonYP, moonZP]);
end


%% =========== Select effective satellites =============
% === satellites are selected base on:
% === 1. The target is assumed at the middle of the surface
% === 2. The signal should not be covered by the Earth
% === 3. The signal should be at least within the side lobe
mainLobeAngle = 23.5;
sideLobeAngle = 90;
targetPosMCMF = llh2MCMF(0, 0, 0);
targetPosECI = mcmf2ECI(targetPosMCMF, moonPosECI, R_MCMF_to_ECI);
targetPosECEF = eci2ecef(utcDate, targetPosECI).';
validSatIndex = [];
mainSatIndex = [];
sideSatIndex = [];
for iSat = 1: numOfSat
    satPosECEF = satPosCurECEF(iSat, :) * 1e3;
    coverFlag = isCovered(satPosECEF, targetPosECEF, earthRadius);
    if coverFlag
        continue;
    end
    lobeFlag = lobeRangeType(satPosECEF, targetPosECEF, mainLobeAngle, sideLobeAngle);
    if lobeFlag ~= 0
        if lobeFlag == 1
            validSatIndex(end + 1) = iSat;
            mainSatIndex(end + 1) = iSat;
        else
            validSatIndex(end + 1) = iSat;
            sideSatIndex(end + 1) = iSat;
        end
    end
end
validSatPos = satPosCurECEF(validSatIndex, :);
mainSatPos = satPosCurECEF(mainSatIndex, :);
sideSatPos = satPosCurECEF(sideSatIndex, :);
numOfValidSat = length(validSatIndex);
% === Plot ===
if enableValidSatPlot
    % === Initialize figure ===
    figure(103);
    hold on;
    grid on;
    xlabel('X (km)', 'FontSize', 18);
    ylabel('Y (km)', 'FontSize', 18);
    zlabel('Z (km)', 'FontSize', 18);
    title('Valid Satellites (ECEF)', 'FontSize', 18);
    view(-58, 13)
    % === Draw Earth ===
    [x, y, z] = sphere(50);
    earthP = surf(x * earthRadius / 1e3, y * earthRadius / 1e3, z * earthRadius / 1e3, ...
        'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', [0 0 1], ...legend
        'DisplayName', 'Earth');
    scale = 1e5;
    earthXP = quiver3(0, 0, 0, scale * x_ECI(1), scale * x_ECI(2), scale * x_ECI(3), ...
        'r', 'LineWidth', 2, 'DisplayName', 'ECI X-axis');
    earthYP = quiver3(0, 0, 0, scale * y_ECI(1), scale * y_ECI(2), scale * y_ECI(3), ...
        'g', 'LineWidth', 2, 'DisplayName', 'ECI Y-axis');
    earthZP = quiver3(0, 0, 0, scale * z_ECI(1), scale * z_ECI(2), scale * z_ECI(3), ...
        'b', 'LineWidth', 2, 'DisplayName', 'ECI Z-axis');
    % === Draw Satellites ===
    % Plot mainlobe sat position
    satMainCP = scatter3(mainSatPos(:, 1), mainSatPos(:, 2), mainSatPos(:, 3), ...
        25, 'o', 'filled', 'MarkerFaceColor', "#77AC30", 'DisplayName', 'Main Lobe');
    % Plot sidelobe sat position
    satSideCP = scatter3(sideSatPos(:, 1), sideSatPos(:, 2), sideSatPos(:, 3), ...
        25, 'o', 'filled', 'MarkerFaceColor', "#0072BD", 'DisplayName', 'Side Lobe');
    % === Draw Moon ===
    moonP = surf((xm * moonRadius + moonPosECEF(1)) / 1e3 , ...
         (ym * moonRadius + moonPosECEF(2)) / 1e3 , ...
         (zm * moonRadius + moonPosECEF(3)) / 1e3, ...
         'FaceAlpha', 1, 'EdgeColor', [0.5 0.5 0.5], 'FaceColor', [0.5 0.5 0.5], ...
         'DisplayName', 'Moon');
    % Scale for visualization
    scale = 5e4;
    % Plot MCMF axes
    moonXP = quiver3(moonPosECEF(1) / 1e3, moonPosECEF(2) / 1e3, moonPosECEF(3) / 1e3, scale * MCMF_X(1), ...
        scale * MCMF_X(2), scale * MCMF_X(3), 'r', 'LineWidth', 2, ...
        'LineStyle', ':', 'DisplayName', 'MCMF X-axis');
    moonYP = quiver3(moonPosECEF(1) / 1e3, moonPosECEF(2) / 1e3, moonPosECEF(3) / 1e3, scale * MCMF_Y(1), ...
        scale * MCMF_Y(2), scale * MCMF_Y(3), 'g', 'LineWidth', 2, ...
        'LineStyle', ':', 'DisplayName', 'MCMF Y-axis');
    moonZP = quiver3(moonPosECEF(1) / 1e3, moonPosECEF(2) / 1e3, moonPosECEF(3) / 1e3, scale * MCMF_Z(1), ...
        scale * MCMF_Z(2), scale * MCMF_Z(3), 'b', 'LineWidth', 2, ...
        'LineStyle', ':', 'DisplayName', 'MCMF Z-axis');

    % === config the figure ===
    axis equal;
    xlim1 = min([0, moonPosECEF(1)]) / 1e3 - 1e5;
    xlim2 = max([0, moonPosECEF(1)]) / 1e3 + 1e5;
    ylim1 = min([0, moonPosECEF(2)]) / 1e3 - 1e5;
    ylim2 = max([0, moonPosECEF(2)]) / 1e3 + 1e5;
    zlim1 = min([0, moonPosECEF(3)]) / 1e3 - 1e5;
    zlim2 = max([0, moonPosECEF(3)]) / 1e3 + 1e5;
    xlim([xlim1, xlim2])
    ylim([ylim1, ylim2])
    zlim([zlim1, zlim2])
    legend([earthP, earthXP, earthYP, earthZP,...
        moonP, moonXP, moonYP, moonZP, satMainCP, satSideCP], 'FontSize', 12);
end



%% =========== Compute GDOP for location selection =============
if enableGDOPPlot
    latCandidates = -90:1:90;
    lonCandidates = -90:1:90; % only the side facing to the Earth
    h = 0;
    GDOP = zeros(length(latCandidates), length(lonCandidates));
    % add extra lunar satellite1
    lunarSatMCMF = llh2MCMF(45, -70, 5200e3);
    lunarSatECI = mcmf2ECI(lunarSatMCMF, moonPosECI, R_MCMF_to_ECI);
    lunarSatECEF(1, :) = eci2ecef(utcDate, lunarSatECI).' / 1e3;

    % lunarSatMCMF = llh2MCMF(45, -45, 5200e3);
    % lunarSatECI = mcmf2ECI(lunarSatMCMF, moonPosECI, R_MCMF_to_ECI);
    % lunarSatECEF(2, :) = eci2ecef(utcDate, lunarSatECI).' / 1e3;

    figure(102)
    plot3(satPosCurECEF(end, 1), satPosCurECEF(end, 2), satPosCurECEF(end, 3), ...
                    '.', 'MarkerSize', 20, 'Color', "#77AC30", 'HandleVisibility', 'off');
    
    for iLat = 1: length(latCandidates)
        for iLon = 1:length(lonCandidates)
            lon = lonCandidates(iLon);
            lat = latCandidates(iLat);
            targetPosMCMF = llh2MCMF(lat, lon, h);
            targetPosECI = mcmf2ECI(targetPosMCMF, moonPosECI, R_MCMF_to_ECI);
            targetPosECEF = eci2ecef(utcDate, targetPosECI).';

            % === select all satellites which are not blocked ===
            % validSatIndex = [];
            % for iSat = 1: numOfSat
            %     satPosECEF = satPosCurECEF(iSat, :) * 1e3;
            %     coverFlag = isCovered(satPosECEF, targetPosECEF, earthRadius);
            %     if coverFlag
            %         continue;
            %     end
            %     validSatIndex(end + 1) = iSat;
            % end
            % validSatPosECEF = satPosCurECEF(validSatIndex,:) * 1e3;
            
            % === select satellites which are valid (beam pattern) ===
            validSatPosECEF = [validSatPos; lunarSatECEF] * 1e3; 

            % === Compute GDOP ===
            GDOP(iLat, iLon) = computeGDOP(validSatPosECEF, targetPosECEF);
        end
    end
    % === visualize GDOP ===
    figure(104)
    % [XPlot, YPlot] = meshgrid(latCandidates, lonCandidates);   % Create a grid for lat and lon
    imagesc(latCandidates, lonCandidates, GDOP);   
    axis xy;  % Correct orientation
    colorbar; % Show color scale
    xlabel('longitude');
    ylabel('latitude');
    title('GDOP on the Moon Surface');
    grid on;
end

%% =========== Additional Functions =============
function E = solveKepler(M, e)
    % SOLVEKEPLER Solves Kepler's equation for eccentric anomaly.
    E = M; % Initial guess
    tol = 1e-8;
    for iter = 1:10
        f = E - e .* sin(E) - M;
        f_prime = 1 - e .* cos(E);
        E_new = E - f ./ f_prime;
        if max(abs(E_new - E)) < tol
            E = E_new;
            break;
        end
        E = E_new;
    end
end

function dydt = equationsOfMotion(~, y, mu)
    % EQUATIONSOFMOTION Computes derivatives of position and velocity.
    r = y(1:3);
    v = y(4:6);
    rNorm = norm(r);
    a = -mu * r / rNorm^3; % Gravitational acceleration
    dydt = [v; a];
end

function posMCMF = llh2MCMF(lat, lon, h)
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

function posECI = mcmf2ECI(posMCMF, l, R)
    % MCMF2ECI Converts MCMF coordinates to ECI coordinates of the Earth
    % xyz - MCMF coordinates
    % l - translation from moon's center to Earth's center
    % R - Rotation Matrix from MCMF to ECI
    tmp = R * posMCMF;
    posECI = l + tmp;
end

function coveredByEarth = isCovered(satPos, targetPos, r)
    % satPos    -- satellite position [x, y, z]
    % targetPos -- target position on the Moon [x, y, z]
    % r         -- radius of the earth

    % Vector from satPos to targetPos
    rho = targetPos - satPos;
       
    % Coefficients for the quadratic equation: At^2 + Bt + C = 0
    A = rho * rho.';
    B = 2 * dot(satPos, rho);
    C = satPos * satPos.' - r^2;
    
    % Discriminant
    delta = B^2 - 4 * A * C;
    
    if delta < 0
        % No intersection
        coveredByEarth = false;
    else
        % Solve for t
        sqrtDelta = sqrt(delta);
        t1 = (-B - sqrtDelta) / (2 * A);
        t2 = (-B + sqrtDelta) / (2 * A);
        
        % Check if either intersection point lies on the segment
        if (t1 >= 0 && t1 <= 1) || (t2 >= 0 && t2 <= 1)
            coveredByEarth = true;  % The line segment intersects the sphere
        else
            coveredByEarth = false; % The infinite line intersects but not the segment
        end
    end
end

function GDOP = computeGDOP(satPos, userPos)
    % satPos: Nx3 matrix of satellite positions [x, y, z]
    % userPos: 1x3 vector of user position [x, y, z]

    N = size(satPos, 1);  % Number of satellites
    A = zeros(N, 4);      % Initialize geometry matrix

    % Construct the geometry matrix A
    for i = 1:N
        dx = satPos(i, 1) - userPos(1);
        dy = satPos(i, 2) - userPos(2);
        dz = satPos(i, 3) - userPos(3);
        
        rho = sqrt(dx^2 + dy^2 + dz^2);  % Distance to satellite
        
        % Fill the geometry matrix A
        A(i, :) = [dx / rho, dy / rho, dz / rho, 1];
    end

    % Compute GDOP
    Q = inv(A' * A);  % Covariance matrix
    GDOP = sqrt(trace(Q));
end

function lobeFlag = lobeRangeType(satPos, userPos, mainLobeAngle, sideLobeAngle)
    % satPos: 1x3 matrix of satellite positions [x, y, z]
    % userPos: 1x3 vector of user position [x, y, z]
    % mainLobeAngle: angle (in degrees) defining the main lobe beamwidth
    % sideLobeAngle: angle (in degrees) defining the side lobe beamwidth
    
    % vector from sat to user on the moon
    a = userPos - satPos;
    % vector of direct main lobe (from sat to the earth)
    b = - satPos;

    % angle of connection between sat and user
    % Compute the angle between vectors a and b
    cosTheta = dot(a, b) / (norm(a) * norm(b));
    theta = acosd(cosTheta); % Convert to degrees

    % Determine which lobe the angle falls into
    if theta <= mainLobeAngle
        lobeFlag = 1;  % Inside main lobe
    elseif theta <= sideLobeAngle
        lobeFlag = 2;  % Between main and side lobe
    else
        lobeFlag = 0;  % Outside side lobe
    end
end

function lobeFlag = LobeRange2(satPos, userPos, mainLobeAngle, sideLobeAngle)
    % satPos: 1x3 matrix of satellite positions [x, y, z]
    % userPos: 1x3 vector of user position [x, y, z]
    % mainLobeAngle: angle (in degrees) defining the main lobe beamwidth
    % sideLobeAngle: angle (in degrees) defining the side lobe beamwidth
end