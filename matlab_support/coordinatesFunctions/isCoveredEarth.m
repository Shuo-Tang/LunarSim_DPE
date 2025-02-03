function coveredByEarth = isCoveredEarth(satPos, targetPos)
    % satPos    -- satellite position [x, y, z]
    % targetPos -- target position on the Moon [x, y, z]
    % r         -- radius of the earth

    r = 6.356752314245179e6;
    % Vector from satPos to targetPos
    rho = targetPos - satPos;
       
    % Coefficients for the quadratic equation: At^2 + Bt + C = 0
    A = rho * rho.';
    B = 2 * dot(satPos, rho);
    C = satPos * satPos.' - r^2;
    
    % Discriminant
    delta = B^2 - 4 * A * C;
    
    if delta < 0
        % No intersection
        coveredByEarth = false;
    else
        % Solve for t
        sqrtDelta = sqrt(delta);
        t1 = (-B - sqrtDelta) / (2 * A);
        t2 = (-B + sqrtDelta) / (2 * A);
        
        % Check if either intersection point lies on the segment
        if (t1 >= 0 && t1 <= 1) || (t2 >= 0 && t2 <= 1)
            coveredByEarth = true;  % The line segment intersects the sphere
        else
            coveredByEarth = false; % The infinite line intersects but not the segment
        end
    end
end