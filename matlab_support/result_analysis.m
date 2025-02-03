clear;
clc;
close all;

%% =========== Load Settings =============
% ==== include necessary packages ====
addpath include             % The software receiver functions
addpath geoFunctions        % Position calculation related functions 
addpath(genpath('goGPS'))
rng(1);
% ==== read settings in yaml ====
settingsFilePath = "../importData/navSettings_1h.yaml";
settings = readyaml(settingsFilePath);
% ==== further settings ====
settings.satRate = 1 / settings.satFrequency;
settings.moonRate = 1 / settings.moonFrequency;
settings.numOfUpdates = settings.satFrequency * settings.timeDuration;
settings.samplesPerMs = settings.samplingFrequency / 1e3;
settings.numOfMs = settings.timeDuration * 1e3;
settings.rlsIter = 10;
settings.numOfSol = settings.timeDuration * settings.solutionFrequency;
settings.solutionRate = 1 / settings.solutionFrequency;
%% =========== Clarify Constants =============
constants = setConstants;
%% =========== Load Satellite Information =============
% ==== read Rinex files ====
ephGPS = [];
ephGlo = [];
ephGal = [];
if settings.enableGPS
    [constellationsGPS] = goGNSS.initConstellation(settings.enableGPS,0,0,0,0,0);
    [ephGPS, ~] = load_RINEX_nav(settings.filePaths.GPS, constellationsGPS,0);
end
if settings.enableGLONASS
    [constellationsGlo] = goGNSS.initConstellation(0,settings.enableGLONASS,0,0,0,0);
    [ephGlo, ~] = load_RINEX_nav(settings.filePaths.GLONASS, constellationsGlo,0);
end
if settings.enableGalileo
    [constellationsGal] = goGNSS.initConstellation(0,0,settings.enableGalileo,0,0,0);
    [ephGal, ~] = load_RINEX_nav(settings.filePaths.Galileo, constellationsGal,0);
end
% ==== compute satellite position ====
% compute GPS time
utcTime0 = datetime(settings.timeEpoch,'TimeZone','UTC');
gpsTime0 = utcDate2gpsSec(utcTime0);
gpsTime = gpsTime0 + (0:1/settings.satFrequency:settings.timeDuration);
% get number of satellites and prns recorded in eph
[~, ~, xs, ~, ~, ~, ~] = satellite_positions(gpsTime(1),...
        ones(32, 36), (1:32), ephGPS, [], [], zeros(32,1), zeros(32,1), 0);
settings.prnGps = find(xs(:,1) ~= 0);
settings.numOfSatGps = length(settings.prnGps);
[~, ~, xs, ~, ~, ~, ~] = satellite_positions(gpsTime(1),...
        ones(24, 36), (1:24), ephGlo, [], [], zeros(24,1), zeros(24,1), 0);
settings.prnGlo = find(xs(:,1) ~= 0);
settings.numOfSatGlo= length(settings.prnGlo);
[~, ~, xs, ~, ~, ~, ~] = satellite_positions(gpsTime(1),...
        ones(36, 36), (1:36), ephGal, [], [], zeros(36,1), zeros(36,1), 0);
settings.prnGal = find(xs(:,1) ~= 0);
settings.numOfSatGal= length(settings.prnGal);

% compute satellite position and velocity
fprintf("Generating satellites position and velocity......\n")
fprintf("GPS......")
[satPosGps, satVelGps] = getSatPosVel('G', settings, gpsTime, ephGPS);
fprintf("GLONASS......")
[satPosGlo, satVelGlo] = getSatPosVel('R', settings, gpsTime, ephGlo);
fprintf("GALILEO......")
[satPosGal, satVelGal] = getSatPosVel('E', settings, gpsTime, ephGal);
fprintf('\n');
fprintf('Done. \n');

%% =========== Load Lunar User Information =============
fprintf("Generating lunar receiver position and velocity......\n")
[moonPosECI, ~] = loadMoonMoving(settings.filePaths.moonPosVel);
moonPosECEF = zeros(size(moonPosECI));
for iMoon = 1:size(moonPosECEF, 1)
    moonPosECEF(iMoon,:) = eci2ecef(utcTime0, moonPosECI(iMoon, :));
end
userLlh = [settings.lunarReceiver.latitude, settings.lunarReceiver.longitude,...
    settings.lunarReceiver.altitude];
[userPosECEF, userVelECEF] = getUserPosVelECEF(userLlh, moonPosECI, utcTime0, settings);
fprintf('Done. \n');


% positioning during 2h
estPos2SP = load('result\estPos2SP_1h.mat');
estPosDPE = load('result\estPosDPE_1h.mat');
estPos2SP = estPos2SP.estPos2SP_2h;
estPosDPE = estPosDPE.estPosDPE_2h;
userPosLlh = [settings.lunarReceiver.latitude,...
    settings.lunarReceiver.longitude, settings.lunarReceiver.altitude];

llhPos2SP = zeros(7200, 3);
enuPos2SP = zeros(7200, 3);
llhPosDPE = zeros(7200, 3);
enuPosDPE = zeros(7200, 3);
% ===== convert to desired coordinates ==== %
for iSol = 1:7200
    eciPos2SP = ecef2eci(utcTime0, estPos2SP(iSol,:));
    llhPos2SP(iSol,:) = eci2moonLlh(eciPos2SP, moonPosECI(iSol,:)')'; 
    enuPos2SP(iSol,:) = moonLlh2enu(llhPos2SP(iSol,:), userPosLlh)';

    eciPosDPE = ecef2eci(utcTime0, estPosDPE(iSol,:));
    llhPosDPE(iSol,:) = eci2moonLlh(eciPosDPE, moonPosECI(iSol,:)')'; 
    enuPosDPE(iSol,:) = moonLlh2enu(llhPosDPE(iSol,:), userPosLlh)';
end
% ==== compute mean estimate and error ==== %
meanPos2SP = zeros(7200, 3);
meanPosDPE = zeros(7200, 3);

for iSol = 2:7200
    meanPos2SP(iSol,:) = mean(enuPos2SP(2:iSol,:), 1);
    meanPosDPE(iSol,:) = mean(enuPosDPE(2:iSol,:), 1);
end

meanError2SP = vecnorm(meanPos2SP, 2, 2);
meanErrorDPE = vecnorm(meanPosDPE, 2, 2);
% ==== plot result ==== %
t = (1:7200) / 60;  % Adjust as needed

% error plot
figure(201);
plot(t, meanError2SP, '-o', 'LineWidth', 1, 'Color', '#d10d0d','DisplayName', '2SP Error');
hold on;
plot(t, meanErrorDPE, '-s', 'LineWidth', 1, 'Color', '#0072bd', 'DisplayName', 'DPE Error');

% Add labels and title
xlabel('Time (minutes)', 'FontSize', 14);
ylabel('Error (meters)', 'FontSize', 14);
title('Positioning Errors Over Time', 'FontSize', 14);
legend;
grid on;

% % ==== compute mean estimate and error ==== %
% meanPos2SP = zeros(7200, 3);
% meanPosDPE = zeros(7200, 3);
% for iSol = 2:7200
%     meanPos2SP(iSol,:) = mean(enuPos2SP(2:iSol,:), 1);
%     meanPosDPE(iSol,:) = mean(enuPosDPE(2:iSol,:), 1);
% end
% meanError2SP = vecnorm(meanPos2SP, 2, 2);
% meanErrorDPE = vecnorm(meanPosDPE, 2, 2);

% pos plot
% Extract East and North components
east2SP = meanPos2SP(2:end,1); % East component of 2SP
north2SP = meanPos2SP(2:end,2); % North component of 2SP

eastDPE = meanPosDPE(2:end,1); % East component of DPE
northDPE = meanPosDPE(2:end,2); % North component of DPE

% Create figure
figure(202);
hold on;

% Normalize indices for color mapping
tColor = linspace(0, 1, length(east2SP)); % Color scale from 0 to 1

% Define custom colormaps
redColors = [linspace(1, 0.4, length(tColor))', linspace(0.8, 0, length(tColor))',...
    linspace(0.8, 0, length(tColor))']; % Light to dark red
blueColors = [linspace(0, 0, length(tColor))', linspace(0.6, 0, length(tColor))',...
    linspace(1, 0.4, length(tColor))']; % Light to dark blue

% Scatter plot with gradually-changing color
scatter(east2SP, north2SP, 50, redColors, 'filled', 'DisplayName', '2SP Solutions');
scatter(eastDPE, northDPE, 50, blueColors, 'filled', 'DisplayName', 'DPE Solutions');

% Add labels and title
xlabel('East (meters)', 'FontSize', 14);
ylabel('North (meters)', 'FontSize', 14);
title('Position Estimate in Lunar ENU Frame', 'FontSize', 14);
legend;
grid on;
axis equal; % Ensures equal scaling

%% show converging result
% finalEst2SP = [meanPos2SP(end,1), meanPos2SP(end,2), meanPos2SP(end,3)];
% finalEstDPE = [meanPosDPE(end,1), meanPosDPE(end,2), meanPosDPE(end,3)];
% 
% fprintf("2D error on E-N plane \n")
% fprintf("2SP:")
% norm(finalEst2SP(1:2))
% fprintf("DPE:")
% norm(finalEstDPE(1:2))
% fprintf("error on height \n")
% fprintf("2SP:")
% norm(finalEst2SP(3))
% fprintf("DPE:")
% norm(finalEstDPE(3))

%% plot with integration time
% err2SP_10_60 = load('result\err2SP_int_time_10_60.mat');
% errDPE_10_60 = load('result\errDPE_int_time_10_60.mat');
% err2SP_70_250 = load('result\err2SP_int_time.mat');
% errDPE_70_250 = load('result\errDPE_int_time.mat');
% err2SP = [err2SP_10_60.err2SP_int_time(:,1:6), err2SP_70_250.err2SP_int_time];
% errDPE = [errDPE_10_60.errDPE_int_time(:,1:6), errDPE_70_250.errDPE_int_time];
% rmse2SP = sqrt(mean(err2SP.^2, 1));
% rmseDPE = sqrt(mean(errDPE.^2, 1));

% err2SP_40_200 = load("result\errPos2SP_ti_40_200.mat");
% errDPE_40_200 = load("result\errPosDPE_ti_40_200.mat");
% err2SP_40_200 = err2SP_40_200.errPos2SP_ti_40_200;
% errDPE_40_200 = errDPE_40_200.errPosDPE_ti_40_200;
% % error plot
% figure(204);
% plot(40:10:200, err2SP_40_200, 'LineWidth', 3, 'Color', '#d10d0d','DisplayName', '2SP Error');
% hold on;
% plot(40:10:200, errDPE_40_200, 'LineWidth', 3, 'Color', '#0072bd', 'DisplayName', 'DPE Error');
% 
% % Add labels and title
% xlabel('Coherent Integration Time (ms)', 'FontSize', 14);
% ylabel('Error (meters)', 'FontSize', 14);
% title('Coherent Integration Time vs Positioning Errors ', 'FontSize', 14);
% legend;
% grid on;