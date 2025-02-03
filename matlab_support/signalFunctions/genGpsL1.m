function CAcode = genGpsL1(prn)
    % Define the G2 shifts for each PRN
    g2s = [5, 6, 7, 8, 17, 18, 139, 140, 141, 251, ...
           252, 254, 255, 256, 257, 258, 469, 470, 471, 472, ...
           473, 474, 509, 512, 513, 514, 515, 516, 859, 860, ...
           861, 862, 145, 175, 52, 21, 237, 235, 886, 657, ...
           634, 762, 355, 1012, 176, 603, 130, 359, 595, 68, ...
           386];

    % Validate PRN
    if prn < 1 || prn > numel(g2s)
        error('Invalid PRN number. Must be between 1 and %d.', numel(g2s));
    end

    % G2 shift for the PRN
    g2shift = g2s(prn);

    % Generate G1 code
    reg1 = -ones(1, 10); % Initial register
    g1 = zeros(1, 1023);
    for i = 1:1023
        g1(i) = reg1(10);
        saveBit = reg1(3) * reg1(10);
        reg1 = [saveBit, reg1(1:9)];
    end

    % Generate G2 code
    reg2 = -ones(1, 10); % Initial register
    g2 = zeros(1, 1023);
    for i = 1:1023
        g2(i) = reg2(10);
        saveBit = prod(reg2([2, 3, 6, 8, 9, 10])); % Feedback taps
        reg2 = [saveBit, reg2(1:9)];
    end

    % Shift G2
    g2Shifted = circshift(g2, g2shift);

    % Generate C/A code
    CAcode = -(g1 .* g2Shifted);
end
