function e1Code = genGalE1(prn)
    % Validate PRN
    if prn < 1 || prn > 36
        error('Invalid PRN number. Must be between 1 and 36.');
    end

    % Initialize G1 and G2 registers
    reg1 = -ones(1, 13);
    reg2 = -ones(1, 13);
    g1 = zeros(1, 4092);
    g2 = zeros(1, 4092);

    % Generate G1 code
    for i = 1:4092
        g1(i) = reg1(13);
        saveBit = prod(reg1([1, 4, 8, 13])); % Feedback taps
        reg1 = [saveBit, reg1(1:12)];
    end

    % Generate G2 code
    for i = 1:4092
        g2(i) = reg2(13);
        saveBit = prod(reg2([1, 2, 9, 13])); % Feedback taps
        reg2 = [saveBit, reg2(1:12)];
    end

    % Combine G1 and G2 to create E1 code
    e1Code = g1 .* circshift(g2, prn);
end
