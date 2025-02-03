function posECI = mcmf2eci(posMCMF, l, R)
    % MCMF2ECI Converts MCMF coordinates to ECI coordinates of the Earth
    % xyz - MCMF coordinates
    % l - translation from moon's center to Earth's center
    % R - Rotation Matrix from MCMF to ECI
    tmp = R * posMCMF;
    posECI = l + tmp;
end