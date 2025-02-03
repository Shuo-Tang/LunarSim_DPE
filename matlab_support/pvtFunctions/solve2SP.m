function [estPos, fdBinGps, fdBinGlo, fdBinGal] =...
                            solve2SP(satGps, satGlo, satGal, ...
                            satPosGpsMs, satPosGloMs, satPosGalMs,...
                            gpsSignal, gloSignal, galSignal, ...
                            localGps, localGlo, localGal, ...
                            fdBinGps, fdBinGlo, fdBinGal,...
                            settings, constants, estPosPrev)

    fracRange = [];
    satPos = [];
    
    if settings.enableGPS
        rGps = zeros(length(satGps), settings.samplesPerMs);
        for iSat = 1:length(satGps)
            prn = settings.prnGps(satGps(iSat));
            [rGps(iSat,:), ~, fdBinGps(prn, :)] = computeCorr('G', gpsSignal, localGps(prn,:), settings, fdBinGps(prn, :));
        end
        fracRangeGps = getFracPseudorange('G', rGps, settings, constants);
        fracRange = [fracRange; fracRangeGps];
        satPos = [satPos; satPosGpsMs];
    end

    if settings.enableGLONASS
        rGlo = zeros(length(satGlo), settings.samplesPerMs);
        for iSat = 1:length(satGlo)
            prn = settings.prnGlo(satGlo(iSat));
            [rGlo(iSat,:), ~, fdBinGlo(prn, :)] = computeCorr('R', gloSignal(prn,:), localGlo(prn,:), settings, fdBinGlo(prn, :));
        end
        fracRangeGlo = getFracPseudorange('R', rGlo, settings, constants);
        fracRange = [fracRange; fracRangeGlo];
        satPos = [satPos; satPosGloMs];
    end

    if settings.enableGalileo
        rGal = zeros(length(satGal), settings.samplesPerMs * 4);
        for iSat = 1:length(satGal)
            prn = settings.prnGal(satGal(iSat));
            [rGal(iSat,:), ~, fdBinGal(prn, :)] = computeCorr('E', galSignal, localGal(prn,:), settings, fdBinGal(prn, :));
        end
        fracRangeGal = getFracPseudorange('E', rGal, settings, constants);
        fracRange = [fracRange; fracRangeGal];
        satPos = [satPos; satPosGalMs];
    end

    estPos = pvtRLS(estPosPrev, satPos, fracRange, settings, constants);
end
