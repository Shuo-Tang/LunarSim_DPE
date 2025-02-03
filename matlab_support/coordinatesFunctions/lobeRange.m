function lobeType = lobeRange(satPos, userPos, mainLobeAngle, sideLobeAngle)
    % satPos    -- satellite position [x, y, z]
    % targetPos -- target position on the Moon [x, y, z]
        

    r = 6.356752314245179e6; % earth radius

    % Vector from satPos to targetPos
    a = userPos - satPos;
       
    % Coefficients for the quadratic equation: At^2 + Bt + C = 0
    A = a * a.';
    B = 2 * dot(satPos, a);
    C = satPos * satPos.' - r^2;
    
    % Discriminant
    delta = B^2 - 4 * A * C;
    
    % LOS covered by Earth or face to the Earth
    if delta > 0
        lobeType = 0;
    % main lobe or side lobe
    else
        b = - satPos;
    
        % angle of connection between sat and user
        % Compute the angle between vectors a and b
        cosTheta = dot(a, b) / (norm(a) * norm(b));
        theta = acosd(cosTheta); % Convert to degrees
    
        % Determine which lobe the angle falls into
        if theta <= mainLobeAngle
            lobeType = 1;  % Inside main lobe
        elseif theta <= sideLobeAngle
            lobeType = 2;  % Between main and side lobe
        else
            lobeType = 0;  % Outside side lobe
        end
    end
    
end

