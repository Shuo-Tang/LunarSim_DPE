function moonPosECI = getMoonPosECI(utcTime)
    % Compute Moon's position a specific epoch defined
    % in UTC time in MATLAB datetime format.
    % Inputs:
    %   utcTime - UTC time in MATLAB datetime format.
    % Output:
    %   moonPosECI   - 1x3 vector of Moon's position in ECI coordinates (meters).
    %
    % include Moon Position package from NASA Jet Propulsion Lab.
    % References:
    % Montenbruck O., Gill E.; Satellite Orbits: Models, Methods and 
    % Applications; Springer Verlag, Heidelberg; Corrected 3rd Printing (2005).
    % Montenbruck O., Pfleger T.; Astronomy on the Personal Computer; Springer 
    % Verlag, Heidelberg; 4th edition (2000).
    % Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill;
    % New York; 4th edition (2013).
    %
    % https://ssd.jpl.nasa.gov/planets/eph_export.html
    %
    % modified by Shuo Tang, NEU to only compute moon's position with 
    % equatorial rectangular coordinates referred to the FK5 equator and equinox 
    % which can be treated as equal to the ECI coordinate system. 
    addpath(genpath('MoonPosition'));
    global const
    SAT_Const
    % get UTC time array
    utcNum = datevec(utcTime);
    % read Earth orientation parameters
    fid = fopen('EOP-All.txt','r');
    %  ----------------------------------------------------------------------------------------------------
    % |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    % |(0h UTC)           "         "          s          s          "        "          "         "     s 
    %  ----------------------------------------------------------------------------------------------------
    while ~feof(fid)
        tline = fgetl(fid);
        k = strfind(tline,'NUM_OBSERVED_POINTS');
        if (k ~= 0)
            numrecsobs = str2num(tline(21:end));
            tline = fgetl(fid);
            for i=1:numrecsobs
                eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
            end
            for i=1:4
                tline = fgetl(fid);
            end
            numrecspred = str2num(tline(22:end));
            tline = fgetl(fid);
            for i=numrecsobs+1:numrecsobs+numrecspred
                eopdata(:,i) = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 1]);
            end
            break
        end
    end
    fclose(fid);

    % compute moon position under ECI (ICRF realized at J2000: Jan 1 2000 12:00:00 TT)
    MJD_UTC = Mjday(utcNum(1), utcNum(2), utcNum(3), utcNum(4), utcNum(5), utcNum(6));
    [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
    [~, ~, ~, TT_UTC, ~] = timediff(UT1_UTC, TAI_UTC);
    MJD_TT = MJD_UTC + TT_UTC/86400;
    % fprintf('very accurate ELP2000-82 lunar coordinates [km]');

    JD_TT = MJD_TT+2400000.5;
    [X,Y,Z] = moonpos(JD_TT);
    moonPosECI = [X,Y,Z].' * 1e3;

end