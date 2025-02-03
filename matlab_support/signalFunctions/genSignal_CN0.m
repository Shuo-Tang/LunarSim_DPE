function [cn0] = genSignal_CN0(constellation, constants, settings, satPos, satVel, userPosECEF, userVelECEF, startMs, endMs)
% GENERATESIGNAL Generate satellite navigation signals for GPS, GLONASS, or Galileo.
%
% This function generates navigation signals for a specified constellation 
% (GPS, GLONASS, or Galileo) using satellite positions, velocities, and 
% user position/velocity data. It computes the modulated signal, Doppler 
% shifts, and fractional delays while checking satellite visibility based on 
% user-defined parameters.
%
% INPUTS:
%   - constellation: (char) Satellite constellation:
%       'G' for GPS
%       'R' for GLONASS
%       'E' for Galileo
%   - constants: (struct) System constants, including:
%       .c             - Speed of light (m/s)
%       .codePeriodGps - Code period for GPS (s)
%       .codePeriodGlo - Code period for GLONASS (s)
%       .codePeriodGal - Code period for Galileo (s)
%       .codeLengthGpsL1 - Code length for GPS L1
%       .codeLengthGloL1 - Code length for GLONASS L1
%       .codeLengthGalE1 - Code length for Galileo E1
%       .TcGps         - Chip period for GPS (s)
%       .TcGlo         - Chip period for GLONASS (s)
%       .TcGal         - Chip period for Galileo (s)
%       .fGpsL1        - Carrier frequency for GPS L1 (Hz)
%       .fGloL1        - Vector of carrier frequencies for GLONASS L1 (Hz)
%       .fGalE1        - Carrier frequency for Galileo E1 (Hz)
%   - settings: (struct) Simulation parameters, including:
%       .numOfUpdates      - Number of updates (integer)
%       .numOfSatGps       - Number of GPS satellites (integer)
%       .numOfSatGlo       - Number of GLONASS satellites (integer)
%       .numOfSatGal       - Number of Galileo satellites (integer)
%       .prnGps            - PRN list for GPS satellites (array)
%       .prnGlo            - PRN list for GLONASS satellites (array)
%       .prnGal            - PRN list for Galileo satellites (array)
%       .samplesPerMs      - Samples per millisecond (integer)
%       .lobeAngleLimit    - Struct with angle limits for main/side lobes
%   - satPos: (matrix) Satellite positions (N x 3), where N is the number of satellites.
%   - satVel: (matrix) Satellite velocities (N x 3), where N is the number of satellites.
%   - userPosECEF: (matrix) User positions in ECEF coordinates (numOfUpdates x 3).
%   - userVelECEF: (matrix) User velocities in ECEF coordinates (numOfUpdates x 3).
%   - fs: (double) Sampling frequency (Hz).
%   - K: (double) Total simulation time in milliseconds.
%   - startMs, endMs: Generate the signal during millisecond [startMs, endMs]
%
% OUTPUTS:
%   - signal: (vector/matrix) Generated signal:
%       For GPS and Galileo: A 1D vector (size (K × fs × 1e-3)).
%       For GLONASS: A 2D matrix (numSat × (K × fs × 1e-3)).
%   - validSat: (cell array) Indices of valid satellites for each update.
%   - local: (matrix) Local replica for future usage

    fs = settings.samplingFrequency;
    K = endMs - startMs + 1;
    % Initialize outputs
    signal = [];
    local = [];
    validSat = cell(K, 1);
    if startMs == endMs
        validSat = [];
    end
    % freqMax = 0;
    % freqMin = 0;
    
    % Determine settings based on the constellation
    switch constellation
        case 'G' % GPS
            codePeriod = constants.codePeriodGps;
            codeLength = constants.codeLengthGpsL1;
            Tc = constants.TcGps;
            fc = constants.fGpsL1;
            numOfSat = settings.numOfSatGps;
            prnList = settings.prnGps;
            genCode = @genCAcode;
            signal = zeros(1, round(K / 1e3 * fs));
            local = zeros(32, round(K / 1e3 * fs));
            cn0 = zeros(32,1);
        case 'R' % GLONASS
            codePeriod = constants.codePeriodGlo;
            codeLength = constants.codeLengthGloL1;
            Tc = constants.TcGlo;
            numOfSat = settings.numOfSatGlo;
            prnList = settings.prnGlo;
            genCode = @genGloL1;
            signal = zeros(24, round(K / 1e3 * fs));
            local = zeros(24, round(K / 1e3 * fs));
            cn0 = zeros(24,1);
        case 'E' % Galileo
            codePeriod = constants.codePeriodGal;
            codeLength = constants.codeLengthGalE1;
            Tc = constants.TcGal;
            fc = constants.fGalE1;
            numOfSat = settings.numOfSatGal;
            prnList = settings.prnGal;
            genCode = @genGalE1;
            signal = zeros(1, round(K / 1e3 * fs));
            local = zeros(36, round(K / 1e3 * fs));
            cn0 = zeros(36,1);
        otherwise
            error('Invalid constellation input. Use "G", "R", or "E".');
    end
    
    % Main loop for signal generation
    for k = startMs:endMs % For each ms
        validSatIndex = [];
        startIdx = 1 + (k - startMs) * settings.samplesPerMs;
        endIdx = (k - startMs  + 1) * settings.samplesPerMs;
        % startIdx = 1 + (k - 1) * settings.samplesPerMs;
        % endIdx = k * settings.samplesPerMs;
        % if startMs == endMs
        %     startIdx = 1;
        %     endIdx = settings.samplesPerMs;
        % end
        
        iUpdate = floor((k - 1) / (1/settings.satFrequency*1e3)) + 1;
        for iSat = 1:numOfSat
            
            satIndex = (iUpdate - 1) * numOfSat + iSat;
            prn = prnList(iSat);
            satPosCurrent = satPos(satIndex, :);
            satVelCurrent = satVel(satIndex, :);
            userPos = userPosECEF(iUpdate, :);
            userVel = userVelECEF(iUpdate, :);

            % Check satellite signal lobe
            lobeFlag = lobeRange(satPosCurrent, userPos, ...
                settings.lobeAngleLimit.mainLobe, settings.lobeAngleLimit.sideLobe);
            range = norm(satPosCurrent - userPos);
            % lobeFlag = isCoveredEarth(satPosCurrent, userPos);
            if lobeFlag == 0
                continue;
            else
                validSatIndex(end + 1) = iSat; %#ok<AGROW>
                CNo = settings.lobeCN0.mainLobe - 20 * log10(range/constants.d) + settings.lobeCN0.mainLobe/10 * randn;% Main lobe
         
                if lobeFlag ~= 1
                    CNo = settings.lobeCN0.sideLobe - 20 * log10(range/constants.d) + settings.lobeCN0.sideLobe/10 * randn; % Side lobe
                end
            end
            
            alpha = sqrt(10^(CNo/10)/settings.samplingFrequency); 

            % record CN0
            cn0(prn) = CNo;
        end 
    end
    
end
