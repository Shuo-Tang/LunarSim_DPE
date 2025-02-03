function [r, cafData] = computeCAF(userPos, satPosMs, satVelMs, userVel, constellation, ...
    constants, settings, cafData, signal, localSignal, dpeFreqBin, satGNSS)
    % Compute CAF contribution for a given GNSS system
    r = 0;
    if constellation == 'G'
        fc = constants.fGpsL1;
        codePeriod = constants.codePeriodGps;
        prnList = settings.prnGps(satGNSS);
    elseif constellation == 'R'
        codePeriod = constants.codePeriodGlo;
        prnList = settings.prnGlo(satGNSS);
    elseif constellation == 'E'
        codePeriod = constants.codePeriodGal;
        fc = constants.fGalE1;
        prnList = settings.prnGal(satGNSS);
    end

    range = vecnorm(satPosMs - userPos, 2, 2);
    unitLos = (satPosMs - userPos) ./ range;
    fracDelay = mod(range / constants.c, codePeriod);
    shift = round(fracDelay * settings.samplingFrequency);
    if constellation ~= 'E'
        shift (shift == settings.samplesPerMs) = 0;
    else
        shift (shift == 4*settings.samplesPerMs) = 0;
    end

    for iSat = 1:length(prnList)
        prn = prnList(iSat);
        if constellation == 'R'
            fc = constants.fGloL1(mod(prn, 14));
        end
        doppler = -((satVelMs(iSat, :) - userVel) * unitLos(iSat, :)') * fc / constants.c;
        [~, fdIdx] = min(abs(dpeFreqBin - doppler));

        
        if isempty(cafData{prn})
            if constellation ~= 'R'
                [~, cafData{prn}, ~] = computeCorr(constellation, signal, localSignal(prn, :), settings, dpeFreqBin);
            else
                [~, cafData{prn}, ~] = computeCorr(constellation, signal(prn, :), localSignal(prn, :), settings, dpeFreqBin);
            end
        end

        r = r + cafData{prn}(fdIdx, shift(iSat) + 1);
    end
end
