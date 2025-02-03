function constants = setConstants()

constants = struct( ...
    'c', 299792458.0, ... % Speed of light (m/s)
    'fGpsL1', 1575.42e6, ... % GPS L1 frequency (Hz)
    'fGalE1', 1575.42e6, ... % Galileo E1 frequency (Hz)
    'fGloL1', [ ... % GLONASS L1 frequencies (Hz) for k = -7 to +6
        1598.0625e6, ... % k = -7
        1598.6250e6, ... % k = -6
        1599.1875e6, ... % k = -5
        1599.7500e6, ... % k = -4
        1600.3125e6, ... % k = -3
        1600.8750e6, ... % k = -2
        1601.4375e6, ... % k = -1
        1602.0000e6, ... % k = 0
        1602.5625e6, ... % k = +1
        1603.1250e6, ... % k = +2
        1603.6875e6, ... % k = +3
        1604.2500e6, ... % k = +4
        1604.8125e6, ... % k = +5
        1605.3750e6  ... % k = +6
    ], ...
    'codePeriodGps', 1.0e-3, ... % GPS code period (s)
    'codePeriodGal', 4.0e-3, ... % Galileo code period (s)
    'codePeriodGlo', 1.0e-3, ... % GLONASS code period (s)
    'TcGps', 1.0 / 1.023e6, ... % GPS chip period (s)
    'TcGal', 4.0 / 4.092e6, ... % Galileo chip period (s)
    'TcGlo', 1.0 / 0.511e6, ... % GLONASS chip period (s)
    'codeLengthGpsL1', 1023, ... % GPS C/A code length
    'codeLengthGalE1', 4092, ... % Galileo E1 code length
    'codeLengthGloL1', 511, ... % GLONASS L1 code length
    'd', 384400e3 ... % reference distance between the moon and the Earth
);



end

