clear;
clc;
close all;

%% =========== Constants / Settings =============
mu = 3.986004418e14; % Earth's gravitational constant (m^3/s^2)
earthRadius = 6.356752314245179e6; %6371e3; % Earth's radius (meters)
moonRadius = 1737e3; % Moon's radius (meters)
OMEGAE_DOT_GPS = 7.2921151467e-5;             % GPS     Angular velocity of the Earth rotation [rad/s]
OMEGAE_DOT_GLO = 7.292115e-5;                 % GLONASS Angular velocity of the Earth rotation [rad/s]
OMEGAE_DOT_GAL = 7.2921151467e-5;             % Galileo Angular velocity of the Earth rotation [rad/s]
OMEGAE_DOT_BDS = 7.292115e-5;                 % BeiDou  Angular velocity of the Earth rotation [rad/s]
CIRCLE_RAD = 2* pi;

settingsFilePath = "../importData/navSettings.yaml";
settings = readyaml(settingsFilePath);
posVelSavePath = "../exportData/moonPosVelMatlab_1Hz_7200s.csv";

%% =========== Computation =============
timeEpoch = settings.timeEpoch;
timeDuration = settings.timeDuration; % Duration in seconds
updateFreq = settings.satFrequency;
updateRate = 1 / updateFreq; 
t0 = datetime(timeEpoch(1), timeEpoch(2), timeEpoch(3), ...
    timeEpoch(4), timeEpoch(5), timeEpoch(6), 0);
numOfUpdates = timeDuration * updateFreq;
tDate = t0 + milliseconds((0: updateRate: timeDuration - updateRate) * 1e3);

moonPosECEF = zeros(numOfUpdates, 3);
% Initialize progress bar
%hWaitBar = waitbar(0, 'Processing epochs...');

parfor iEpoch = 1:numOfUpdates
    iEpoch
    % fprintf("Porcessing: %d / %d /n", iEpoch, numOfUpdates)
    % Compute Moon's position
    moonPosECI = getMoonPosECI(tDate(iEpoch));
    moonPosECEF(iEpoch, :) = eci2ecef(tDate(iEpoch), moonPosECI);
    % Update progress bar
    % waitbar(iEpoch / numOfUpdates, hWaitBar, sprintf('Processing epoch %d of %d...', iEpoch, numOfUpdates));
end
% Estimate Moon's velocity
moonVelECEF = diff(moonPosECEF, 1, 1);
moonVelECEF = [moonVelECEF(1,:);moonVelECEF];

moonPosVelECEF = [moonPosECEF, moonVelECEF];

% Close the progress bar
%close(hWaitBar);
disp('Processing complete.');

% save result
writematrix(moonPosVelECEF, posVelSavePath);