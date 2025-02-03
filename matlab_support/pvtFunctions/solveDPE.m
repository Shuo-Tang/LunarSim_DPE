function [estPos, cafGps, cafGlo, cafGal] = solveDPE(satPosGpsMs, satVelGpsMs, satPosGloMs, satVelGloMs, ...
                  satPosGalMs, satVelGalMs, userPos, userVel, ...
                  cafGps, cafGlo, cafGal, gpsSignal, gloSignal, galSignal, ...
                  localGps, localGlo, localGal, satGps, satGlo, satGal, ...
                  dpeFreqBin, settings, constants)

    % extract settings
    dMax = settings.ARS.dmax;
    dMin = settings.ARS.dmin;
    decay = settings.ARS.decay;
    nIter = settings.ARS.nIter;
    gammaEst = zeros(settings.ARS.nIter, 3);
    J = zeros(settings.ARS.nIter, 1);
    % set initial point
    d = dMax; 
    gamma = [userPos(1:2) + 1000*(2*rand(1,2)-1), userPos(3)];
    % gamma = userPos + d*(2*rand(1,3)-1);
    gammaEst(1,:) = gamma;
    r = 0;
    if settings.enableGPS
        [rGps, cafGps] = computeCAF(gamma, satPosGpsMs, satVelGpsMs, ...
            userVel, 'G', constants, settings,...
            cafGps, gpsSignal, localGps, dpeFreqBin, satGps);
        r = r + rGps;
    end
    if settings.enableGLONASS
        [rGlo, cafGlo] = computeCAF(gamma, satPosGloMs, satVelGloMs, ...
            userVel, 'R', constants, settings,...
            cafGlo, gloSignal, localGlo, dpeFreqBin, satGlo);
        r = r + rGlo;
    end
    if settings.enableGalileo
        [rGal, cafGal] = computeCAF(gamma, satPosGalMs, satVelGalMs, ...
            userVel, 'E', constants, settings,...
            cafGal, galSignal, localGal, dpeFreqBin, satGal);
        r = r + rGal;
    end
    % compute initial cost
    cost = r;
    J(1) = cost;
    % ars iteration
    for it = 1:nIter-1
        % draw a random movement
        % gamma = gammaEst(it,:) + d*(2*rand(1,3)-1);
        gamma = [gammaEst(it,1:2) + d*(2*rand(1,2)-1), userPos(3)];
        % compute cost
        r = 0;
        if settings.enableGPS
            [rGps, ~] = computeCAF(gamma, satPosGpsMs, satVelGpsMs, ...
                userVel, 'G', constants, settings,...
                cafGps, gpsSignal, localGps, dpeFreqBin, satGps);
            r = r + rGps;
        end
        if settings.enableGLONASS
            [rGlo, ~] = computeCAF(gamma, satPosGloMs, satVelGloMs, ...
                userVel, 'R', constants, settings,...
                cafGlo, gloSignal, localGlo, dpeFreqBin, satGlo);
            r = r + rGlo;
        end
        if settings.enableGalileo
            [rGal, ~] = computeCAF(gamma, satPosGalMs, satVelGalMs, ...
                userVel, 'E', constants, settings,...
                cafGal, galSignal, localGal, dpeFreqBin, satGal);
            r = r + rGal;
        end
        % select or discard point
        if r > cost
            gammaEst(it+1,:) = gamma;
            cost = r;
            d = dMax;
        else
            gammaEst(it+1,:) = gammaEst(it,:);
            d = d * decay;
        end
        J(it+1) = cost;
        if d < dMin
            d = dMax;
        end
    end

    estPos = gammaEst(end, :);

    % ==== plot search path ====
    % figure(201); hold on;
    % path = plot3(gammaEst(:,1) - userPos(1), gammaEst(:,2) - userPos(2), J,...
    %     'LineWidth',3, 'Color',"#A2142F");
    % legend(path, "ARS Path", 'FontSize',14)

end