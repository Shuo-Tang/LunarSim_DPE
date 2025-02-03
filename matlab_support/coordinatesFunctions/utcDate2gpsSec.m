function gpsTime = utcDate2gpsSec(utcTime)
%UTCDATE2GPSSEC corrects an array of UTC dates(in any matlab format) into
%GPS time in seconds counting from Jan 6 1980 00:00:00
%   Shuo Jan 2025

gpsDate0 = datetime([1980, 1, 6, 0, 0, 0],'TimeZone','UTC');
gpsTimeSerial = utc2gps(utcTime);
gpsDate = datetime(gpsTimeSerial, "ConvertFrom", "datenum",'TimeZone','UTC');
gpsTime = seconds(gpsDate - gpsDate0);
end

