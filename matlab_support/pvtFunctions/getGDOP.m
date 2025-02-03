function GDOP = getGDOP(satPos, userPos)
    % satPos: Nx3 matrix of satellite positions [x, y, z]
    % userPos: 1x3 vector of user position [x, y, z]

    N = size(satPos, 1);  % Number of satellites
    A = zeros(N, 4);      % Initialize geometry matrix

    % Construct the geometry matrix A
    for i = 1:N
        dx = satPos(i, 1) - userPos(1);
        dy = satPos(i, 2) - userPos(2);
        dz = satPos(i, 3) - userPos(3);
        
        rho = sqrt(dx^2 + dy^2 + dz^2);  % Distance to satellite
        
        % Fill the geometry matrix A
        A(i, :) = [dx / rho, dy / rho, dz / rho, 1];
    end

    % Compute GDOP
    Q = inv(A' * A);  % Covariance matrix
    GDOP = sqrt(trace(Q));
end