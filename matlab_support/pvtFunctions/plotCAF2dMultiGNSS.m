function plotCAF2dMultiGNSS( ...
    satPosGpsMs, satVelGpsMs, satPosGloMs, satVelGloMs, ...
    satPosGalMs, satVelGalMs, userPos, userVel, ...
    settings, constants,...
    dpeFreqBin, cafGps, gpsSignal, localGps, ...
    cafGlo, gloSignal, localGlo, cafGal, galSignal, localGal, ...
    satGps, satGlo, satGal)

    % Function to plot 2D Cross-Ambiguity Function (CAF) for GPS, GLONASS, and Galileo
    
    % ==== Define search grid ==== %
    density = 2;
    numOfPoints = 500;

    xCan = userPos(1) - numOfPoints * density : density : userPos(1) + numOfPoints * density;
    yCan = userPos(2) - numOfPoints * density : density : userPos(2) + numOfPoints * density;
    r = zeros(length(xCan), length(yCan));

    % ==== Compute CAF over the search grid ==== %
    for ix = 1:length(xCan)
        for iy = 1:length(yCan)
            userPosCan = [xCan(ix), yCan(iy), userPos(3)];
            
            % GPS Processing
            if settings.enableGPS
                [rGps, cafGps] = computeCAF(userPosCan, satPosGpsMs, satVelGpsMs, ...
                    userVel, 'G', constants, settings,...
                    cafGps, gpsSignal, localGps, dpeFreqBin, satGps);
                r(ix, iy) = r(ix, iy) + rGps;
            end

            % GLONASS Processing
            if settings.enableGLONASS
                [rGlo, cafGlo] = computeCAF(userPosCan, satPosGloMs, satVelGloMs, ...
                    userVel, 'R', constants, settings,...
                    cafGlo, gloSignal, localGlo, dpeFreqBin, satGlo);
                r(ix, iy) = r(ix, iy) + rGlo;
            end

            % Galileo Processing
            if settings.enableGalileo
                [rGal, cafGal] = computeCAF(userPosCan, satPosGalMs, satVelGalMs, ...
                    userVel, 'E', constants, settings,...
                    cafGal, galSignal, localGal, dpeFreqBin, satGal);
                r(ix, iy) = r(ix, iy) + rGal;
            end
        end
    end

    % ==== Plot the results ==== %
    [xPlot, yPlot] = meshgrid(xCan - userPos(1), yCan - userPos(2));
    figure(201);
    mesh(xPlot, yPlot, r');
    title('Multi-constellation CAF', 'FontSize',14);
    xlabel('X Offset (m)', 'FontSize',14);
    ylabel('Y Offset (m)', 'FontSize',14);
    zlabel('CAF Value', 'FontSize',14);
end