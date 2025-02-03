function g1 = genGloL1(prn)
    % Define the polynomial G(X) = 1 + X^5 + X^9
    % The corresponding feedback taps for a 9-bit LFSR
    polyTaps = [9, 5];  % Feedback from X^9 and X^5

    % Number of chips per period
    numChips = 511;

    % Initialize LFSR with a nonzero seed (typically all ones)
    LFSR = ones(1, 9);

    % Preallocate memory for the output ranging code
    g1 = zeros(1, numChips);

    % Generate the ranging code
    for i = 1:numChips
        % Output the first bit as the next chip
        g1(i) = LFSR(end);

        % Compute new bit using feedback taps (XOR of selected bits)
        newBit = mod(sum(LFSR(polyTaps)), 2);

        % Shift register to the right and insert new bit at the start
        LFSR = [newBit, LFSR(1:end-1)];
    end

    % Replace 0s with -1s
    g1(g1 == 0) = -1;
end


