function [estFracRange] = getFracPseudorange(constellation, r, settings, constants)
    %% Load configuration
    fs = settings.samplingFrequency;
    c = constants.c;

    %% Memory allocation

    %% Estimate time delays
    [~, maxPos] = max(r, [], 2);
    maxPos = maxPos - 1;
    EstFracDelay = maxPos / fs;
    if constellation == 'E'
        EstFracDelay = mod(EstFracDelay, 1e-3); 
    end
    estFracRange = EstFracDelay * c;
end
