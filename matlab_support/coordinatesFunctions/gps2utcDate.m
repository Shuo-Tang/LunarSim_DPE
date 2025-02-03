function utcDate = gps2utcDate(gpsTime)
    % GPSSECONDSTOUTC Converts GPS seconds to UTC datetime.
    % leap seconds
    leapEpoch = [46828800, 78364801, 109900802, 173059203, 252028804, 315187205,...
                346723206, 393984007, 425520008, 457056009, 504489610, 551750411,...
                599184012, 820108813, 914803214, 1025136015, 1119744016, 1167264017];
    leapSec = 0;
    i = 1;
    while i <= 18 && gpsTime > leapEpoch(i)
        leapSec = leapSec + 1;
        i = i + 1;
    end
    % compute UTC date
    gpsDate = datetime(1980, 1, 6, 0, 0, 0, 'TimeZone', 'UTC');
    utcDate = gpsDate + seconds(gpsTime - leapSec);
end

