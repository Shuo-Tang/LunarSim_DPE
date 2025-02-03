%% This script is the base file of all the other scripts with differnet configurations

clear;
clc;
close all;

%% =========== Load Settings =============
% ==== include necessary packages ====
addpath include                 % The software receiver functions
addpath geoFunctions            % Position calculation related functions 
addpath coordinatesFunctions    % Functions for frame conversion and geometry relationship
addpath pvtFunctions            % Functions for solving 2SP and DPE
addpath signalFunctions         % Functions for signal generation
addpath(genpath('goGPS'))
rng(1);
% ==== read settings in yaml ====
settingsFilePath = "../importData/navSettings.yaml";
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
[moonPosECI, ~] = loadMoonMoving(settings.filePaths.moonPosVel);
userLlh = [settings.lunarReceiver.latitude, settings.lunarReceiver.longitude,...
    settings.lunarReceiver.altitude];
fprintf("Generating lunar receiver position and velocity......\n")
[userPosECEF, userVelECEF] = getUserPosVelECEF(userLlh, moonPosECI, utcTime0, settings);
fprintf('Done. \n');

% =========== Generate Signal and solve PVT =============
fprintf("Solving PVT for the first Epoch......\n")
%% ==== initialzie signal in first integration time ====
if settings.enableGPS
    [gpsSignal, satGps, localGps] = genSignal('G', constants, settings,...
        satPosGps, satVelGps, userPosECEF, userVelECEF, 1, settings.coherentIntegrationTime);
end
if settings.enableGLONASS
    [gloSignal, satGlo, localGlo] = genSignal('R', constants, settings,...
        satPosGlo, satVelGlo, userPosECEF, userVelECEF, 1, settings.coherentIntegrationTime);
end
if settings.enableGalileo
    [galSignal, satGal, localGal] = genSignal('E', constants, settings,...
        satPosGal, satVelGal, userPosECEF, userVelECEF, 1, settings.coherentIntegrationTime);
end
% ==== initialzie positioning solution ====
initialFreqBin = -2e4:500:2e4;
fdBinGps = zeros(32, 8);
fdBinGlo = zeros(24, 8);
fdBinGal = zeros(36, 8);
fracRangeGps = [];
fracRangeGlo = [];
fracRangeGal = [];
if settings.enableGPS
    rGps = zeros(length(satGps{1}), settings.samplesPerMs);
    % get pseudorange measurement
    for iSat = 1:length(satGps{1})
        prn = settings.prnGps(satGps{1}(iSat));
        localSignalGps = localGps(prn,:);
        [rGps(iSat,:), cafGps, fdBinGps(prn, :)] = computeCorr('G', gpsSignal, localSignalGps, settings, initialFreqBin);
    end
    fracRangeGps = getFracPseudorange('G', rGps, settings, constants);
    % compute PVT
    % estPVT = pvtRLS(userPosECEF(1,:), satPosGps(satGps{1},:), fracRangeGps, settings, constants);
end
if settings.enableGLONASS
    rGlo = zeros(length(satGlo{1}), settings.samplesPerMs);
    % get pseudorange measurement
    for iSat = 1:length(satGlo{1})
        prn = settings.prnGlo(satGlo{1}(iSat));
%         localSignalGlo = localGlo(prn,:);
        [rGlo(iSat,:), cafGlo, fdBinGlo(prn, :)] = computeCorr('R', gloSignal(prn,:), localGlo(prn,:), settings, initialFreqBin);
    end
    fracRangeGlo = getFracPseudorange('R', rGlo, settings, constants);
    % compute PVT
    % estPVT = pvtRLS(userPosECEF(1,:), satPosGlo(satGlo{1},:), fracRangeGlo, settings, constants);
end
if settings.enableGalileo
    rGal = zeros(length(satGal{1}), settings.samplesPerMs*4);
    % get pseudorange measurement
    for iSat = 1:length(satGal{1})
        prn = settings.prnGal(satGal{1}(iSat));
        localSignalGal = localGal(prn,:);
        [rGal(iSat,:), cafGal, fdBinGal(prn, :)] = computeCorr('E', galSignal, localSignalGal, settings, initialFreqBin);
    end
    fracRangeGal = getFracPseudorange('E', rGal, settings, constants);
    % compute PVT
    % estPVT = pvtRLS(userPosECEF(1,:), satPosGal(satGal{1},:), fracRangeGal, settings, constants);
end
fracRange = [fracRangeGps; fracRangeGlo; fracRangeGal];
satPos = [satPosGps(satGps{1},:); satPosGlo(satGlo{1},:); satPosGal(satGal{1},:)];
estPVT0 = pvtRLS(userPosECEF(1,:), satPos, fracRange, settings, constants);
fprintf('Done. \n');


%% ==== Generating signal in real time and solve PVT ====
% fprintf("Solving PVT in real time......\n")
% h = waitbar(0, 'Processing...'); % Create waitbar
% totalIterations = settings.numOfMs - settings.coherentIntegrationTime;
% estPos = zeros(settings.numOfMs - settings.coherentIntegrationTime + 1, 3);
% estPos(1,:) = estPVT0; 
% for iMs = 2:(settings.numOfMs - settings.coherentIntegrationTime + 1)
%     % Update waitbar
%     waitbar((iMs - 1) / totalIterations, h, sprintf('Processing %d of %d', iMs - 1, totalIterations));
% 
%     iUpdate = floor((iMs - 1) / (1/settings.satFrequency*1e3)) + 1;
%     fracRangeGps = [];
%     fracRangeGlo = [];
%     fracRangeGal = [];
%     if settings.enableGPS
%         [gpsSignal1ms, satGps1ms, ~] = genSignal('G', constants, settings,...
%             satPosGps, satVelGps, userPosECEF, userVelECEF, iMs, iMs);
%         gpsSignal = [gpsSignal(settings.samplesPerMs + 1:end), gpsSignal1ms];
%         rGps = zeros(length(satGps1ms), settings.samplesPerMs);
%         % get pseudorange measurement
%         for iSat = 1:length(satGps1ms)
%             prn = settings.prnGps(satGps1ms(iSat));
%             localSignalGps = localGps(prn,:);
%             [rGps(iSat,:), cafGps, fdBinGps(prn,:)] = computeCorr('G', gpsSignal, localSignalGps, settings, fdBinGps(prn,:));
%         end
%         fracRangeGps = getFracPseudorange('G', rGps, settings, constants);
%         satIndexGps = (iUpdate - 1) * settings.numOfSatGps + satGps1ms;
%         % compute PVT
%         % estPVT = pvtRLS(estPVT, satPosGps(satIndexGps,:), fracRangeGps, settings, constants);
%     end
%     if settings.enableGLONASS
%         [gloSignal1ms, satGlo1ms, ~] = genSignal('R', constants, settings,...
%             satPosGlo, satVelGlo, userPosECEF, userVelECEF, iMs, iMs);
%         rGlo = zeros(length(satGlo1ms), settings.samplesPerMs);
%         for iSat = 1:length(satGlo1ms)
%             prn = settings.prnGlo(satGlo1ms(iSat));
%             gloSignal(prn, :) = [gloSignal(prn, settings.samplesPerMs + 1:end),...
%                 gloSignal1ms(prn, :)];
%             % localSignalGlo = localGlo(prn,:);
%             [rGlo(iSat,:), cafGlo, fdBinGlo(prn,:)] = computeCorr('R', gloSignal(prn,:),localGlo(prn,:),settings, fdBinGlo(prn,:));
%         end
%         fracRangeGlo = getFracPseudorange('R', rGlo, settings, constants);
%         satIndexGlo = (iUpdate - 1) * settings.numOfSatGlo + satGlo1ms;
%         % estPVT = pvtRLS(estPVT, satPosGlo(satIndexGlo,:), fracRangeGlo, settings, constants);
%     end
%     if settings.enableGalileo
%         [galSignal1ms, satGal1ms, ~] = genSignal('E', constants, settings,...
%             satPosGal, satVelGal, userPosECEF, userVelECEF, iMs, iMs);
%         galSignal = [galSignal(settings.samplesPerMs + 1:end), galSignal1ms];
%         rGal = zeros(length(satGal1ms), 4*settings.samplesPerMs);
%         for iSat = 1:length(satGal1ms)
%             prn = settings.prnGal(satGal1ms(iSat));
%             localSignalGal = localGal(prn,:);
%             [rGal(iSat,:), cafGal, fdBinGal(prn,:)] = computeCorr('E', galSignal, localSignalGal, settings,fdBinGal(prn,:));
%         end
%         fracRangeGal = getFracPseudorange('E', rGal, settings, constants);
%         satIndexGal = (iUpdate - 1) * settings.numOfSatGal + satGal1ms;
%         % estPVT = pvtRLS(estPVT, satPosGal(satIndexGal,:), fracRangeGal, settings, constants);
%     end
%     fracRange = [fracRangeGps; fracRangeGlo; fracRangeGal];
%     satPos = [satPosGps(satIndexGps, :); satPosGlo(satIndexGlo, :); satPosGal(satIndexGal, :);];
%     estPos(iMs,:) = pvtRLS(estPos(iMs-1,:), satPos, fracRange, settings, constants);
% end
% 
% close(h);
% fprintf('Done. \n');



%% ==== Generating signal based on solution rate and solve PVT ====
% ==== Generate signal ==== %
fprintf("Solving PVT based on solution rate......\n")
rng(1);
estPos2SP = zeros(settings.numOfSol, 3);
estPos2SP(1,:) = estPVT0;
errPos2SP = zeros(settings.numOfSol, 1);
errPos2SP(1) = norm(estPVT0 - userPosECEF(1,:));
estPosDPE = zeros(settings.numOfSol, 3);
estPosDPE(1,:) = estPVT0;
errPosDPE = zeros(settings.numOfSol, 1);
errPosDPE(1) = norm(estPVT0 - userPosECEF(1,:));
cafGps = cell(32, 1);
cafGlo = cell(24, 1);
cafGal = cell(36, 1);
h = waitbar(0, 'Solving PVT'); % Create waitbar
for iSol = 2: settings.numOfSol
    % Update waitbar
    waitbar((iSol - 1) / settings.numOfSol, h, sprintf('Solving PVT %d of %d', iSol - 1, settings.numOfSol));
    msStart = (iSol - 1)*settings.solutionRate*1e3 + 1;
    msEnd = msStart + settings.coherentIntegrationTime - 1;
    iUpdate = floor((msStart - 1) / (1/settings.satFrequency*1e3)) + 1;
    if settings.enableGPS
        [gpsSignal, satGps, localGps] = genSignal('G', constants, settings,...
            satPosGps, satVelGps, userPosECEF, userVelECEF, msStart, msEnd);
    end
    if settings.enableGLONASS
        [gloSignal, satGlo, localGlo] = genSignal('R', constants, settings,...
            satPosGlo, satVelGlo, userPosECEF, userVelECEF, msStart, msEnd);
    end
    if settings.enableGalileo
        [galSignal, satGal, localGal] = genSignal('E', constants, settings,...
            satPosGal, satVelGal, userPosECEF, userVelECEF, msStart, msEnd);
    end

    % ==== extract satelliet and user information ==== %
    satIndexGps = (iUpdate - 1) * settings.numOfSatGps + satGps{1};
    satPosGpsMs = satPosGps(satIndexGps,:);
    satVelGpsMs = satVelGps(satIndexGps,:);
    satIndexGlo = (iUpdate - 1) * settings.numOfSatGlo + satGlo{1};
    satPosGloMs = satPosGlo(satIndexGlo,:);
    satVelGloMs = satVelGlo(satIndexGlo,:);
    satIndexGal = (iUpdate - 1) * settings.numOfSatGal + satGal{1};
    satPosGalMs = satPosGal(satIndexGal,:);
    satVelGalMs = satVelGal(satIndexGal,:);
    userPos = userPosECEF(iUpdate, :);
    userVel = userVelECEF(iUpdate, :);

    %% 2SP
    if settings.enable2SP   
        [estPos2SP(iSol,:), fdBinGps, fdBinGlo, fdBinGal]...
                            = solve2SP(satGps{1}, satGlo{1}, satGal{1}, ...
                            satPosGpsMs, satPosGloMs, satPosGalMs,...
                            gpsSignal, gloSignal, galSignal, ...
                            localGps, localGlo, localGal, ...
                            fdBinGps, fdBinGlo, fdBinGal,...
                            settings, constants, estPos2SP(iSol-1,:));
        errPos2SP(iSol) = norm(estPos2SP(iSol,:) - userPos);
    end

    %% DPE
    if settings.enableDPE
        % ==== compute DPE correlegram ==== %
        dpeFreqBin = -2e4:500:2e4;
        numFreqBins = length(dpeFreqBin);
        
        % ==== plot caf ==== %
        plotCAF2dMultiGNSS( ...
        satPosGpsMs, satVelGpsMs, satPosGloMs, satVelGloMs, satPosGalMs, satVelGalMs, ...
        userPos, userVel, ...
        settings, constants,...
        dpeFreqBin, cafGps, gpsSignal, localGps, ...
        cafGlo, gloSignal, localGlo, cafGal, galSignal, localGal, ...
        satGps{1}, satGlo{1}, satGal{1})
    
        % ==== DPE ARS ==== %
        [estPosDPE(iSol, :), cafGps, cafGlo, cafGal] =...
                  solveDPE(satPosGpsMs, satVelGpsMs, satPosGloMs, satVelGloMs, ...
                  satPosGalMs, satVelGalMs, userPos, userVel, ...
                  cafGps, cafGlo, cafGal, gpsSignal, gloSignal, galSignal, ...
                  localGps, localGlo, localGal, satGps{1}, satGlo{1}, satGal{1}, ...
                  dpeFreqBin, settings, constants);
        errPosDPE(iSol) = norm(estPosDPE(iSol,:) - userPos);
       
    end

 end
close(h);
fprintf('Done. \n');
