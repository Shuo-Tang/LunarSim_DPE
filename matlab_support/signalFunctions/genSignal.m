function [signal, validSat, local] = genSignal(constellation, constants, settings, satPos, satVel, userPosECEF, userVelECEF, startMs, endMs)
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
            
        case 'R' % GLONASS
            codePeriod = constants.codePeriodGlo;
            codeLength = constants.codeLengthGloL1;
            Tc = constants.TcGlo;
            numOfSat = settings.numOfSatGlo;
            prnList = settings.prnGlo;
            genCode = @genGloL1;
            signal = zeros(24, round(K / 1e3 * fs));
            local = zeros(24, round(K / 1e3 * fs));
            
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
                CNo = settings.lobeCN0.mainLobe - 20 * log10(range/constants.d);% Main lobe
                if lobeFlag ~= 1
                    CNo = settings.lobeCN0.sideLobe - 20 * log10(range/constants.d); % Side lobe
                end
            end
            alpha = sqrt(10^(CNo/10)/settings.samplingFrequency); 

            % Compute signal parameters
            % range = norm(satPosCurrent - userPos);
            unitLos = (satPosCurrent - userPos) / range;
            fracDelay = mod(range / constants.c, codePeriod);
            % ncoIndex = -fracDelay / Tc;
            shift = round(fracDelay * fs);
            if constellation == 'R'
                fc = constants.fGloL1(mod(prn, 14));
            end

            doppler = -((satVelCurrent - userVel) * unitLos') * fc / constants.c;
            % if doppler > freqMax
            %     freqMax = doppler;
            % end
            % if doppler < freqMin
            %     freqMin = doppler;
            % end

            % Generate code
            code = genCode(prn);
            if constellation == 'G'|| constellation == 'R'
                ii = 1:settings.samplesPerMs;
            else
                ii = 1:4*settings.samplesPerMs;
                iii = 1:settings.samplesPerMs;
            end
                        
            if constellation == 'G'
                % Accumulate the signal
                % signal(startIdx:endIdx) = signal(startIdx:endIdx) + ...
                %     alpha * code(mod(round(ncoIndex + (ii - 1) / fs / Tc), codeLength) + 1);% .* ...
                    %exp(1i * 2 * pi * doppler * (ii - 1) / fs);
                signal(startIdx:endIdx) = signal(startIdx:endIdx) + ...
                    alpha * circshift(code(mod(round((ii - 1) / fs / Tc), codeLength) + 1),shift).* ...
                    exp(1i * 2 * pi * doppler * (ii - 1) / fs);
            elseif constellation == 'R'
                % seperate the signal in prn channels
                % signal(prn, startIdx:endIdx) = ...
                %     alpha * code(mod(round(ncoIndex + (ii - 1) / fs / Tc), codeLength) + 1) .* ...
                %     exp(1i * 2 * pi * doppler * (ii - 1) / fs);
                signal(prn, startIdx:endIdx) = ...
                    alpha * circshift(code(mod(round((ii - 1) / fs / Tc), codeLength) + 1),shift) .* ...
                    exp(1i * 2 * pi * doppler * (ii - 1) / fs);
            elseif constellation == 'E'  
                % code1ms = code(1 + mod(k-1, 4)*1023:(mod(k-1, 4)+1)*1023);
                code4msShift = circshift(code(mod(round((ii - 1) / fs / Tc), codeLength) + 1),shift);
                signal(startIdx:endIdx) = signal(startIdx:endIdx) + ...
                alpha * code4msShift(1 + mod(k-1, 4)*settings.samplesPerMs:(mod(k-1, 4)+1)*settings.samplesPerMs).* ...
                    exp(1i * 2 * pi * doppler * (iii - 1) / fs);
            end
           
            % generate local here to avoid repeating operations
            if constellation == 'G' || constellation == 'R'
                local(prn, startIdx:endIdx) = ...
                    code(mod(round((ii - 1) / fs / Tc), codeLength) + 1);%.* ...
                    % exp(1i * 2 * pi * doppler * (ii - 1) / fs);
            elseif constellation == 'E'
                code4ms = code(mod(round((ii - 1) / fs / Tc), codeLength) + 1);
                local(prn, startIdx:endIdx) = ...
                    code4ms(1 + mod(k-1, 4)*settings.samplesPerMs:(mod(k-1, 4)+1)*settings.samplesPerMs);%.* ...
                    % exp(1i * 2 * pi * doppler * (iii - 1) / fs);
            end

        end
        
        % record valid satelites
        if startMs == endMs
            validSat = validSatIndex;
        else
            validSat{k - startMs + 1} = validSatIndex;
        end
    end
    % add noise
    if constellation ~= 'R'
        noise = (sqrt(1/2)*randn(1,K*settings.samplesPerMs) + 1i*sqrt(1/2)*randn(1,K*settings.samplesPerMs));
        signal = signal + noise;
    else
        for iSat = validSatIndex
            noise = (sqrt(1/2)*randn(1,K*settings.samplesPerMs) + 1i*sqrt(1/2)*randn(1,K*settings.samplesPerMs));
            signal(prnList(iSat),:) = signal(prnList(iSat),:) + noise;
        end
    end
    
end
